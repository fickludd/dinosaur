package se.lth.immun

import Dinosaur._
import DinoUtil._

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet

import se.lth.immun.mzml.ghost.GhostSpectrum

import akka.actor._

case class BuildFromSpectra(n:Int)
case class Ms1Count(n:Int)
case class BuiltHills(hills:Seq[Hill])
case class AddSpectrum(scanIndex:Int, ms1Index:Int, cs:CentSpectrum)

object HillBatch {
	def apply(ms1Range:Ms1Range, hills:Seq[Hill]):HillBatch = {
		val (borderL, temp) = hills.partition(_.scanIndex.head <= ms1Range.start+1)
		val (mid, borderR) = temp.partition(_.lastIndex < ms1Range.end - 2)
		new HillBatch(borderL, mid, borderR)
	}
}

class HillBatch(
		val borderL:Seq[Hill], 
		val mid:Seq[Hill], 
		val borderR:Seq[Hill]
) {
	def toSeq = borderL ++ mid ++ borderR
	
	def merge(hb:HillBatch)(implicit params:DinosaurParams) = {
		val t0 = System.currentTimeMillis
		val ret = 
			new HillBatch(
				borderL, 
				mid ++ sew(borderR.sortBy(_.centerMz.last), hb.borderL.sortBy(_.centerMz.head), params).filter(
					_.scanIndex.length >= params.adv.hillMinLength
				) ++ hb.mid,
				hb.borderR
			)
		if (params.verbose)
			println("%8d %8d %8d %8d %8dms".format(mid.length, borderR.length, hb.borderL.length, hb.mid.length, System.currentTimeMillis - t0))
		ret
	}
		
	def sew(rhs1:Seq[Hill], rhs2:Seq[Hill], params:DinosaurParams):Seq[Hill] = {
		var i1 = 0
		var i2 = 0
		val res = new ArrayBuffer[Hill]
		while (i1 < rhs1.length && i2 < rhs2.length) {
			val rh1 = rhs1(i1)
			val rh2 = rhs2(i2)
			val mz1 = rh1.centerMz.last
			val mz2 = rh2.centerMz.head
			if (rh2.scanIndex.head - rh1.lastIndex <= (params.adv.hillMaxMissing.value+1) && within(mz1, mz2, params.adv.hillPPM.value)) {
				res += rh1.append(rh2)
				i1 += 1
				i2 += 1
			} else if (mz1 < mz2) {
				res += rh1
				i1 += 1
			} else {
				res += rh2
				i2 += 1
			}
		}
		while (i1 < rhs1.length) {
			res += rhs1(i1)
			i1 += 1
		}
		while (i2 < rhs2.length) {
			res += rhs2(i2)
			i2 += 1
		}
		res
	}
}



class HillBuilderActorParallel(val params:DinosaurParams) extends Actor {
	
	implicit val pd = params
	
	var numSpec = -1
	var ms1Count = -1
	var customer:ActorRef = _
	val specMap = new HashMap[Int, Array[CentSpectrum]]
	val specIndexes = new HashMap[Int, Int]
	val toProcess = new Queue[(Int, Array[CentSpectrum])]
	val beingProcessed = new HashSet[Int]
	val done = new HashMap[Ms1Range, Seq[Hill]]
	
	val readProfiler = new ReadProfiler(params.verbose)
	
	def receive = {
		case BuildFromSpectra(n:Int) =>
			customer = sender
			numSpec = n
			val str = readProfiler.initiate("reading file (nSpec=%d) ...".format(numSpec))
			print(str)
		
			
		case AddSpectrum(scanIndex, ms1Index, c) =>
			val start = (ms1Index / params.adv.hillBatchSize) * params.adv.hillBatchSize
			if (start == ms1Index)
				specIndexes(start) = scanIndex
			if (!specMap.contains(start))
				specMap += start -> new Array[CentSpectrum](params.adv.hillBatchSize)
			specMap(start)(ms1Index - start) = c
			checkAndQueue(start)
			processIfFree
			
				
		case BuiltRange(range, rawHills) =>
			done(range) = rawHills
			beingProcessed -= range.start
			processIfFree
			
			print(readProfiler.onRange(specIndexes(range.start), range, rawHills))
			
			
			if (rawHillsDone) {
				if (params.verbose)
					println("sewing started")
				val sewnHills = sewAndRefineHills.toSeq
				
				customer ! BuiltHills(sewnHills)
				
				specMap.clear
				done.clear
				context.stop(self)
			}
			
			
		case Ms1Count(n) =>
			ms1Count = n
			/*
			val start = ((ms1Count-1) / params.adv.hillBatchSize) * params.adv.hillBatchSize
			if (specMap.nonEmpty) {
				if (specMap.keys.toSeq.contains(start)) {
					specMap(start) = specMap(start).take(ms1Count - start)
					checkAndQueue(start)
					processIfFree
				}
			}*/
			
			specMap.size match {
				case 0 => {} // we're good to go, just wait for final processing
				case 1 => 
					val finalBatch = specMap.keys.head
					specMap(finalBatch) = specMap(finalBatch).filter(_ != null) 
					queueForProcess(finalBatch)
					processIfFree
				case n =>
					throw new Exception("Too many hill-batches left after mzML reading complete. Expected 0 or 1, got"+n)
			}
	}
	
	def checkAndQueue(start:Int) = {
		val specs = specMap(start)
		if (specs.forall(_ != null))
			queueForProcess(start)
	}
	
	def queueForProcess(start:Int) = {
		toProcess += start -> specMap(start)
		specMap -= start
	}
	
	def processIfFree = {
		while (toProcess.nonEmpty && beingProcessed.size < params.concurrency) {
			val (start, specs) = toProcess.dequeue
			process(start, specs)
		}
	}
	
	def process(start:Int, specs:Array[CentSpectrum]) = {
		beingProcessed += start
		val a = context.actorOf(Props(new RawHillActor(params)))
		a ! Build(Ms1Range(start, start + specs.length), specs)
	}
	
	def rawHillsDone = 
		ms1Count >= 0 && specMap.isEmpty && toProcess.isEmpty && beingProcessed.isEmpty
		
	def sewAndRefineHills = {
		val rawHills = done.toSeq.sortBy(_._1.start)
		if (params.verbose)
			println("  left n     l-hem   r-hem   right n      time")
		
		val sewn = merge(rawHills)
		if (params.verbose)
			println("sewn raw hills, n="+sewn.length)
		val t0 = System.currentTimeMillis
		
		val runt = Runtime.getRuntime
		if (params.verbose) {
			println("refining hills")
			println("      id     smooth    decompose   lenFilter  massCalc     time     heap")
		}
		
		(for ((batch,i) <- sewn.grouped(params.freqHillRefine).zipWithIndex) yield {
			
			val tt0 = System.currentTimeMillis
			batch.foreach(_.smooth)
			val tt1 = System.currentTimeMillis
			val decomposed = batch.flatMap(_.decompose)
			val tt2 = System.currentTimeMillis
			val filtered = decomposed.filter(_.length >= params.adv.hillMinLength)
			val tt3 = System.currentTimeMillis
			val averaged = filtered.map(_.calcAverages)
			val tt4 = System.currentTimeMillis
			
			if (params.verbose)
				println("%8d %8dms %8dms %8dms %8dms %8ds %8dMb".format(
						i*params.freqHillRefine, 
						tt1 - tt0,
						tt2 - tt1,
						tt3 - tt2,
						tt4 - tt3, 
						(System.currentTimeMillis - t0) / 1000,
						runt.totalMemory / 1000000
					))
			
			averaged
		}).flatten
	}
		
	def merge(hills:Seq[(Ms1Range, Seq[Hill])]):Seq[Hill] = {
		val hillBatches = hills.sortBy(_._1.start).map(t => HillBatch(t._1, t._2))
		var batch = hillBatches.head
		for (hb <- hillBatches.tail)
			batch = batch.merge(hb)
		batch.toSeq
	}
	
	
	class ReadProfiler(active:Boolean) {
		
		var t0 = 0L
		var tTemp = 0L
		
		def onRange(specIndex:Int, range:Ms1Range, rawHills:Seq[Hill]) = {
			if (active) {
				val nScans = range.end - range.start
				val runt = Runtime.getRuntime
				
				val row = "%8d %8d %3d %8d | %8.1f %8.1f | %8dms %7ds %8dMb\n".format(
						specIndex,
						range.start,
						range.end - range.start,
						rawHills.map(_.centerMz.length).sum / nScans,
						
						rawHills.length.toDouble / nScans,
						rawHills.map(_.centerMz.length).sum.toDouble / rawHills.length,
						
						(System.currentTimeMillis - tTemp),
						(System.currentTimeMillis - t0)/1000,
						runt.totalMemory / 1000000
					)
				tTemp = System.currentTimeMillis
				row
			} else ""
		}
		
		def initiate(batch:String) = {
			
			t0 = System.currentTimeMillis
			tTemp = t0
			if (active)	batch + "\n    spec      ms1 ms1n   hLen/scan  nHills   hLenAvg        dt       t       heap\n"
			else ""
		}
	}
}



case class Ms1Range(start:Int, end:Int)
case class Build(range:Ms1Range, cspecs:Array[CentSpectrum])
case class BuiltRange(range:Ms1Range, rawHills:Seq[Hill])
class RawHillActor(val params:DinosaurParams) extends Actor {
	
	var activeHills:Seq[Hill] = Nil
	val completedHills = new ArrayBuffer[Hill] 
	def finalHills(range:Ms1Range) = {
		val all = completedHills ++ activeHills
		all.filter(h => h.lastIndex == range.end-1 || h.scanIndex.head == range.start || h.scanIndex.length >= params.adv.hillMinLength)
	} 
	
	case class ToProcess(ms1Index:Int, cspec:CentSpectrum)
	
	val specQueue = new Queue[ToProcess]
	
	
	
	def receive = {
		case Build(range, specs) =>
			
			for (i <- 0 until specs.length)
				processSpec(range.start+i, specs(i))
			
			sender ! BuiltRange(range, finalHills(range))
			terminate
	}
	
	
	
	def terminate = {
		context.stop(self)
	}
	
	
	
	def processSpec(ms1Index:Int, spec:CentSpectrum) = {
		val mergedHills = merge(activeHills, spec.cpeaks, ms1Index)
		val (active, done) = mergedHills.partition(h => 
								ms1Index - h.lastIndex <= params.adv.hillMaxMissing
							)
		completedHills ++= done
		activeHills = active
	}
	
	
	case class HillPeakDiff(h:Hill, peak:CentroidPeak, mzDiff:Double)
	def merge(hills:Seq[Hill], cpeaks:Seq[CentroidPeak], ms1Index:Int):Seq[Hill] = {
		var ih = 0
		var ic = 0
		val res = new ArrayBuffer[Hill]
		var cluster = new HashSet[HillPeakDiff]
		def resolve = {
			val resolved = new ArrayBuffer[Hill]
			while (cluster.nonEmpty) {
				val min = cluster.minBy(_.mzDiff)
				cluster -= min
				resolved += min.h.push(min.peak, ms1Index)
				cluster = cluster.filter(hpd => hpd.h != min.h)
			}
			res ++= resolved.sortBy(_.mzGuess)
		}
		
		while (ih < hills.length && ic < cpeaks.length) {
			val h = hills(ih)
			val cp = cpeaks(ic)
			val hmz = h.mzGuess
			if (params.closeToDebug(hmz)) 
				{ val k = 1 }
			
			if (within(hmz, cp.mz, params.adv.hillPPM)) {
				cluster += HillPeakDiff(h, cp, math.abs(hmz - cp.mz))
				ih += 1
				ic += 1
			} else {
				resolve
				if (cp.mz < hmz) {
					res += Hill(cp, ms1Index)(params)
					ic += 1
				} else {
					res += h
					ih += 1
				}
			}
		}
		while (ih < hills.length) {
			res += hills(ih)
			ih += 1
		}
		while (ic < cpeaks.length) {
			res += Hill(cpeaks(ic), ms1Index)(params)
			ic += 1
		}
		res
	}
	
}