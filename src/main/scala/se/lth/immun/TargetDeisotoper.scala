package se.lth.immun

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue
import scala.collection.mutable.HashSet

import akka.actor._

import se.lth.immun.chem.Peptide
import se.lth.immun.chem.Constants
import se.lth.immun.chem.IsotopeDistribution

class TargetDeisotoper(val params:DinosaurParams) extends Actor with Timeable {

	import Cluster._
	
	implicit val p = params
	
	val toProcess 			= new Queue[(Int, Seq[Target])]
	val beingProcessed 		= new HashSet[Int]
	val completedPatterns	= new ArrayBuffer[Seq[IsotopePattern]]
	
	var hillsByMz:Array[Hill] = _
	var hillsMzs:Array[Double] = _
	var specTime:Seq[Double] = _
	var customer:ActorRef = _
	
	val targetBatchSize = 1000
	
	def receive = {
		case TargetDeisotope(hs, targets, st) => 
			timer.start
			
			val hills = hs.toArray
			specTime = st
			customer = sender
			
			println("=== IN TARGETED MODE - BEWARE!!! ===")
		  hillsByMz = hills.sortBy(_.total.centerMz)
		  hillsMzs = hillsByMz.map(h => h.total.centerMz)
		
		  if (params.verbose)
			  println("deisotoping based on targets...")
			
			/*
			 * Using a list for the hills is extremely slow for many hills
			 * an array is asympotically faster
			 */
      val batchTargetList = targets.grouped(targetBatchSize).toList
      for (i <- 0 until batchTargetList.length)
			  toProcess += i -> batchTargetList(i)
			
			if (params.verbose)
  			println("created " + batchTargetList.length + " target deisotoping batches")
			processIfFree
		
		case Deisotoped(batchId, isotopes) =>
		  completedPatterns += isotopes
		  beingProcessed -= batchId
		  if (params.verbose)
  		  println("batch deisotoping finished " + batchId)	
		  
		  if (toProcess.isEmpty && beingProcessed.isEmpty) {
				if (params.verbose)
				  println("deisotoping finished")
				
				customer ! TargetDeisotopingComplete(Nil, completedPatterns.flatten)
				context.stop(self)
			} else
			  processIfFree
						
	}
	
	
	def processIfFree = {
		while (toProcess.nonEmpty && beingProcessed.size < params.concurrency) {
			val (batchId, batchTargets) = toProcess.dequeue
			beingProcessed += batchId
			val a = context.actorOf(Props(new TargetBatchDeisotoper(params)))
  		a ! BatchDeisotope(batchId, hillsByMz, batchTargets, hillsMzs, specTime)
		}
	}
}

class TargetBatchDeisotoper(val params:DinosaurParams) extends Actor {
  
  import Cluster._
  
  implicit val p = params
  
  def receive = {
		
		case BatchDeisotope(batchId, hillsByMz, batchTargets, hillsMzs, specTime) =>
		  val (_, isotopes) = deisotope(hillsByMz, batchTargets, hillsMzs, specTime)
		  
		  sender ! Deisotoped(batchId, isotopes)
		  context.stop(self)
  }
  
	def deisotope(hillsByMz:Array[Hill],
			targets:Seq[Target],
			hillsMzs:Array[Double],
			specTime:Seq[Double]			
	):(Seq[Seq[Cluster.Edge]], Seq[IsotopePattern]) = {
		
		val targetPatterns = 
			for (t <- targets) yield {
				val monoisoHills = closeHills(hillsByMz, t, hillsMzs, specTime)
				if (monoisoHills.nonEmpty) {
					val patterns = getPatterns(monoisoHills, hillsByMz, t).map(ip =>
						IsotopePattern(ip.inds.map(hillsByMz), ip.offset, ip.mostAbundNbr - ip.offset, ip.z, ip.averagineCorr))
					
					(t, patterns)
				} else
					(t, Nil)
			}
		
		(Nil, targetPatterns.flatMap(_._2))
	}
	
	
	def closeHills(hills:Array[Hill], t:Target, hillsMzs:Array[Double], specTime:Seq[Double]):Seq[Int] = {
	  val (hillsStartIndx, hillsEndIndx) = DinoUtil.getMinxMaxIndx(hillsMzs, t.mz, t.mzDiff)
		val inds = new ArrayBuffer[Int]
		for (i <- hillsStartIndx until hillsEndIndx) {
			val h = hills(i)
			val hApexRt = h.accurateApexRt(specTime)
			if (hApexRt > t.rtStart && hApexRt < t.rtEnd)
				inds += i
		}
		inds
	}
	
	
	
	def getPatterns(
			seeds:Seq[Int], 
			hills:Array[Hill],
			t:Target
	):Seq[IsotopePatternInds] = {
		if (seeds.isEmpty) return Nil
		val patterns = 
			for {
				seed <- seeds
				isopatInds <- extendSeed(seed, Nil, hills, t.z)
			} yield (isopatInds)
		if (patterns.isEmpty) Nil
		else 
			patterns.filter(ip => ip.inds.head == ip.seed)
	}
	
	
	def extendSeed(seed:Int, inds:Seq[Int], hills:Array[Hill], z:Int):Option[IsotopePatternInds] = {
		
		val seedTot = hills(seed).total
			
		def extend2(dir:Int):Seq[Int] = {
			var ii = seed+dir
			var nIso = 1
			val sHill = hills(seed)
			val sTot = sHill.total
			var m = sTot.centerMz + dir * nIso * DinoUtil.ISOTOPE_PATTERN_DIFF / z
			var seedErr2 = sTot.centerMzError * sTot.centerMzError
			var isoMissing = false
			val isos = new ArrayBuffer[Int]
			val alts = new ArrayBuffer[Int]
			
			
			def evalAlts = {
				if (alts.nonEmpty) {
					val corrs = alts.map(a => (a,params.adv.deisoCorrCalc(hills(seed), hills(a))))
					val maxCorr = corrs.maxBy(_._2)
					alts.clear
					if (maxCorr._2 >= params.adv.deisoCorr) {
						isos += maxCorr._1
						nIso += 1
						m = seedTot.centerMz + dir * nIso * DinoUtil.ISOTOPE_PATTERN_DIFF / z
					} else 
						isoMissing = true
				} else
					isoMissing = true
			}
			
			while (!isoMissing && ii < hills.length && ii >= 0) {
				val iHill = hills(ii)
			  val iTot = iHill.total
			  val mDiff = iTot.centerMz - m
			  val massErrorSq = params.adv.deisoSigmas * params.adv.deisoSigmas * (
					  seedErr2 + iTot.centerMzError * iTot.centerMzError)
			  val err2 = DinoUtil.SULPHUR_SHIFT * DinoUtil.SULPHUR_SHIFT / (z*z) + massErrorSq
			  if (mDiff * mDiff <= err2) {
			    if (overlap(sHill, iHill) > 0) {
  				  alts += ii
  				}
				  ii += dir
			  } else if (dir*mDiff > 0) {
				  evalAlts
			  } else {
				  ii += dir
			  }
			}
			evalAlts
			isos
		}
					
		val upMatches = extend2(1)
		//val downMatches = extend2(-1)
		val downMatches = Nil
		
		val result = downMatches.reverse ++ (seed +: upMatches)
		val resSeedInd = downMatches.length
		val resultProfile = result.map(hills(_).apex.intensity)
		
		val minima = DinoUtil.localMinima(resultProfile, params.adv.deisoValleyFactor)
		val oneMaxResult = 
			if (minima.nonEmpty) {
				val lower = minima.filter(_ < resSeedInd).lastOption.getOrElse(0)
				val upper = minima.filter(_ > resSeedInd).headOption.getOrElse(result.length)
				result.slice(lower, upper + 1)
			} else result
		
		val cleanResult =
			if (z * seedTot.centerMz < 1000) {
				val apex = oneMaxResult.maxBy(hills(_).total.intensity)
				oneMaxResult.drop(oneMaxResult.indexOf(apex))
			} else oneMaxResult
		
		val cleanProfile = cleanResult.map(hills(_).total.intensity)
		val avgIsotopeDistr = Peptide.averagine((seedTot.centerMz - Constants.PROTON_WEIGHT)*z).getIsotopeDistribution()
		
		val alignment = params.adv.deisoAveragineStrategy(cleanProfile, avgIsotopeDistr, params)
		
		alignment.map(a => IsotopePatternInds(
				cleanResult,//.drop(math.max(0, -a.offset)), 
				0, //math.max(0, a.offset), 
				avgIsotopeDistr.intensities.indexOf(avgIsotopeDistr.intensities.max), 
				z, 
				a.corr,
				seed
			))
	}
	
	
	
	def overlapping(seed:Int, hills:Array[Hill]):Seq[Int] =
		for {
			i <- 0 until hills.length
			if overlap(hills(seed), hills(i)) > 0
		} yield i
	
	
	
	
	def overlap(h1:Hill, h2:Hill) =
		math.min(h1.scanIndex.last, h2.scanIndex.last) - math.max(h1.scanIndex.head, h2.scanIndex.head) 
}

case class BatchDeisotope(batchId:Int, hillsByMz:Array[Hill], batchTargets:Seq[Target], hillsMzs:Array[Double],specTime:Seq[Double])
case class Deisotoped(batchId:Int, isotopes:Seq[IsotopePattern])
case class TargetDeisotope(hills:Seq[Hill], targets:Seq[Target], specTime:Seq[Double])
case class TargetDeisotopingComplete(clusters:Seq[Seq[Cluster.Edge]], patterns:Seq[IsotopePattern])
