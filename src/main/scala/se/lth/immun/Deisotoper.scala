package se.lth.immun

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue
import scala.collection.mutable.HashSet
import scala.collection.mutable.HashMap

import akka.actor._

import Deisotoper._
import DinoUtil._

import se.lth.immun.chem.Peptide
import se.lth.immun.chem.Constants
import se.lth.immun.chem.IsotopeDistribution

object Deisotoper {
	
	
	class EdgeProfiler(val repFreq:Int, active:Boolean) {
		var pairCount = 0
		var corrCount = 0
		var okCount = 0
		var totalCount = 0
		var profileLength = 0
		var totalTime = new Timer
		var batchTime = new Timer
		
		def reportAndReset(i:Int) = {
			totalCount += okCount
			println("%8d %8d %8.3f | %6.1f %6.1f | %7ds %6d %8dms".format(
					i, 
					pairCount / repFreq, 
					corrCount.toDouble / repFreq,
					profileLength.toDouble / repFreq,
					okCount.toDouble / repFreq,
					totalTime.check / 1000, 
					totalCount,
					batchTime.click
				))
			pairCount = 0
			corrCount = 0
			profileLength = 0
			okCount = 0
		}
		
		def initiate(batch:String) = {
			if (active) {
				println(batch)
				println("base hill  edge/h   r2Comp/h   hLen   ok/h      tot t  tot ok     time")
			}
			totalTime.start
			batchTime.start
			totalCount = 0
		}
		
		def onHill(i:Int) = 
			if (active && i % repFreq == 0)
				reportAndReset(i)
	
	}
}


case class Deisotope(hills:Seq[Hill], specTime:Seq[Double])
case class DeisotopingComplete(clusters:Seq[Seq[Cluster.Edge]], patterns:Seq[IsotopePattern])

class Deisotoper(val params:DinosaurParams) extends Actor with Timeable {
	
	import Cluster._
	
	implicit val pd = params
	
	val toProcess 			= new Queue[(Int, Seq[Seq[Edge]])]
	val beingProcessed 		= new HashSet[Int]
	val completedPatterns	= new ArrayBuffer[Seq[IsotopePattern]]
	val clusters 			= new ArrayBuffer[Seq[Edge]]
	
	var hills:Seq[Hill] = _
	var customer:ActorRef = _ 
		
	
	
	
	/*
	 * hills need to be sorted by maxRT / scanIndex
	 */
	def receive = {
		case Deisotope(hs, specTime) => 
			timer.start
			
			hills = hs
			customer = sender
			
			/*
			 * Using a list for the hills is extremely slow for many hills
			 * an array is asympotically faster
			 */
			val ahills = hills.toArray
			val profiler = new EdgeProfiler(params.freqFindEdge, params.verbose)
			
			params.findDebugHill(hills, specTime)
			
			profiler.initiate("find egdes")
			val edges = findEdges(ahills, profiler, params.debugHillInd)
			
			println("edges assembled")
			params.deisoEdgeTime = timer.click
			
			def ripEdges(i:Int):Seq[Edge] = {
				val q = Queue[Int](i)
				val res = new ArrayBuffer[Edge]
				while (q.nonEmpty) {
					val j = q.dequeue
					val es = edges(j)
					if (es != null) {
						edges(j) = null
						res ++= es
						for (e <- es)
							if (edges(e.j) != null)
								q += e.j
					}
				}
				res
				/*
				 * This initial recursive impl causes stack overflow for some files
				if (edges(i) != null) {
					val es = edges(i)
					edges(i) = null
					es ++ es.flatMap(e => ripEdges(e.j))
				}
				else Nil //List(Edge(i, i, 0, 0))
				*/
			}
			
			
			for (i <- 0 until edges.length) {
				val cluster = ripEdges(i)
				if (cluster.nonEmpty)
					clusters += cluster
			}
			
			params.findDebugCluster(clusters)
			
			params.deisoRipTime = timer.click
			println("clusters ripped, n="+clusters.length)
			for ((clusts, i) <- clusters.grouped(params.freqClusterDeconvolve).zipWithIndex) 
				toProcess += i -> clusts
			
			if (params.verbose)
				println(" cluster   isotopes   comp time      time")
			processIfFree
						
			
			
		case Deconvolved(id, isotopes, t) =>
			completedPatterns += isotopes
			beingProcessed -= id
			
			if (params.verbose)
				println("%8d %8d | %8dms %8ds".format(
					id * params.freqClusterDeconvolve,
					isotopes.length,
					t,
					(System.currentTimeMillis - timer.t0)/1000
				))
			
			if (toProcess.isEmpty && beingProcessed.isEmpty) {
				println("clusters deconvolved")
				val t3 = System.currentTimeMillis
				
				params.deisoDeconvolveTime = timer.click
				
				customer ! DeisotopingComplete(clusters, completedPatterns.flatten)
				context.stop(self)
			} else 
				processIfFree
	}
	
	
	
	
	def processIfFree = {
		while (toProcess.nonEmpty && beingProcessed.size < params.concurrency) {
			val (id, clusters) = toProcess.dequeue
			beingProcessed += id
			val a = context.actorOf(Props(new Deconvolver(params)))
			a ! Deconvolve(id, clusters, hills)
		}
	}
	
	
	
	def findEdges(hills:Seq[Hill], profiler:EdgeProfiler, stopHillIndex:Option[Int]):Array[List[Edge]] = {
		
		val edges = new Array[List[Edge]](hills.length)
		def addEdge(i:Int, j:Int, z:Int, corr:Double) = {
			if (edges(i) == null) 
				edges(i) = Nil
			edges(i) = Edge(i, j, z, corr) :: edges(i)
		}
		
		for (i <- 1 until hills.length) {
			profiler.onHill(i)
			val hi = hills(i)
			val minScanIndex = hi.scanIndex.head
			val minInd = hills.lastIndexWhere(_.scanIndex.last < minScanIndex, i)+1
			profiler.pairCount += i - minInd
			for (j <- minInd until i) {
				val hj = hills(j)
				val mDiff = math.abs(hi.total.centerMz - hj.total.centerMz)
				if (stopHillIndex.map(x => x==i || x == j).getOrElse(false))
					{ val k = 1 }
				if (mDiff <= C13C12_DIFF + 0.1) {
					
					val mErrPart = hi.total.centerMzError * hi.total.centerMzError + hj.total.centerMzError * hj.total.centerMzError
					val mError = 5 * math.sqrt(mErrPart)
					
					if (mDiff <= C13C12_DIFF + mError) {
						fitsDiff(mDiff, mError) match {
							case Some(z) =>
								val corr = params.adv.deisoCorrCalc(hi, hj)
								profiler.corrCount += 1
								profiler.profileLength += math.max(hi.scanIndex.last, hj.scanIndex.last) - math.min(hi.scanIndex.head, hj.scanIndex.head)
								if (corr > params.adv.deisoCorr) {
									profiler.okCount += 1
									addEdge(i, j, z, corr)
									addEdge(j, i, z, corr)
								}
							case None => {}
						}
					}
				}
			}
		}
		edges
	}
	
	def fitsDiff(mDiff:Double, mError:Double):Option[Int] = {
		for (z <- params.minCharge.value to params.maxCharge) {
			val error = math.sqrt(SULPHUR_SHIFT * SULPHUR_SHIFT / (z*z) + mError*mError)
			if (math.abs(mDiff - ISOTOPE_PATTERN_DIFF / z) <= error)
				return Some(z)
		}
		None
	}
}




case class Deconvolve(id:Int, edgeLists:Seq[Seq[Cluster.Edge]], hills:Seq[Hill])
case class Deconvolved(id:Int, patterns:Seq[IsotopePattern], t:Long)

class Deconvolver(val params:DinosaurParams) extends Actor {
	
	implicit val pd = params
	import Cluster._
	
	def deconvolve(edgeLists:Seq[Seq[Edge]], hills:Seq[Hill]) = 
		edgeLists.flatMap(edgeList => Cluster(edgeList, hills).deconvolve())
	
	
	def receive = {
		case Deconvolve(id, edgeLists, hills) =>
			val t0 = System.currentTimeMillis
			val patterns = deconvolve(edgeLists, hills)
			sender ! Deconvolved(id, patterns, System.currentTimeMillis - t0)
			context.stop(self)
		
	}
}