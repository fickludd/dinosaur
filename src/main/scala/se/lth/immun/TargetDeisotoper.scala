package se.lth.immun

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue
import scala.collection.mutable.HashSet

import akka.actor._

import se.lth.immun.chem.Peptide
import se.lth.immun.chem.Constants
import se.lth.immun.chem.IsotopeDistribution

object TargetDeisotoper {
  def filterUnique(isotopes:Seq[IsotopePattern], specTime:Seq[Double]) = 
    isotopes.foldRight(List.empty[IsotopePattern]){ 
      case (a, b) => 
        if (!b.isEmpty && b(0).hills.head.total.centerMz == a.hills.head.total.centerMz && b(0).z == a.z && b(0).apexHill.accurateApexRt(specTime) == a.apexHill.accurateApexRt(specTime)) 
          b 
        else 
          a :: b 
    }
}

class TargetDeisotoper(val params:DinosaurParams) extends Actor with Timeable {

	import Cluster._
	
	implicit val p = params
	
	val toProcess 			= new Queue[(Int, Seq[Target])]
	val beingProcessed 		= new HashSet[Int]
	val completedPatterns	= new ArrayBuffer[Seq[IsotopePattern]]
	
	var hillsByMz:Array[Hill] = _
	var hillsMzs:Array[Double] = _
	var hillsRts:Array[Double] = _
	var hillsScanIndices:Array[(Int, Int)] = _
	var specTime:Seq[Double] = _
	var customer:ActorRef = _
	
	val targetBatchSize = 1000
	
	def receive = {
		case TargetDeisotope(hs, targets, st) => 
			timer.start
			
			/*
			 * Using a list for the hills is extremely slow for many hills
			 * an array is asympotically faster
			 */
			val hills = hs.toArray
			specTime = st
			customer = sender
			
			println("=== IN TARGETED MODE ===")
		  hillsByMz = hills.sortBy(_.total.centerMz)
		  hillsMzs = hillsByMz.map(h => h.total.centerMz)
		  hillsRts = hillsByMz.map(h => h.accurateApexRt(specTime))
		  hillsScanIndices = hillsByMz.map(h => (h.scanIndex.head, h.scanIndex.last))
		  
		  if (params.verbose)
			  println("deisotoping based on targets...")

      val batchTargetList = targets.sortBy(_.mz).grouped(targetBatchSize).toList
      for ((batch, i) <- batchTargetList.zipWithIndex)
			  toProcess += i -> batch
			
			if (params.verbose)
  			println("created " + batchTargetList.length + " target deisotoping batches")
			processIfFree
		
		case Deisotoped(batchId, isotopes) =>      
		  completedPatterns += isotopes
		  beingProcessed -= batchId
		  if (params.verbose)
  		  println("batch deisotoping finished " + batchId + ", found " + isotopes.length + " unique isotope patterns")	
		  
		  if (toProcess.isEmpty && beingProcessed.isEmpty) {
				if (params.verbose)
				  println("deisotoping finished")
				
				val isotopesAll = completedPatterns.flatten.sortBy(ip => (ip.hills.head.total.centerMz, ip.z, ip.apexHill.accurateApexRt(specTime)))
				if (params.verbose)
				  println("sorting finished, " + isotopesAll.length + " isotope patterns found")
				
				val uniqueIsotopes = TargetDeisotoper.filterUnique(isotopesAll, specTime)
				
        if (params.verbose)
				  println("filtered unique isotope patterns, " + uniqueIsotopes.length + " isotope patterns remaining")
				
				customer ! TargetDeisotopingComplete(Nil, uniqueIsotopes.toSeq)
				context.stop(self)
			} else
			  processIfFree
						
	}
	
	
	def processIfFree = {
		while (toProcess.nonEmpty && beingProcessed.size < params.concurrency) {
			val (batchId, batchTargets) = toProcess.dequeue
			beingProcessed += batchId
			val a = context.actorOf(Props(new TargetBatchDeisotoper(params)))
  		a ! BatchDeisotope(batchId, hillsByMz, batchTargets, hillsMzs, hillsRts, hillsScanIndices, specTime)
		}
	}
}

class TargetBatchDeisotoper(val params:DinosaurParams) extends Actor {
  
  import Cluster._
  
  implicit val p = params
  
  def receive = {
		
		case BatchDeisotope(batchId, hillsByMz, batchTargets, hillsMzs, hillsRts, hillsScanIndices, specTime) =>
		  val (_, isotopes) = deisotope(hillsByMz, batchTargets, hillsMzs, hillsRts, hillsScanIndices)
		  	  
		  val uniqueIsotopes = TargetDeisotoper.filterUnique(isotopes.sortBy(ip => (ip.hills.head.total.centerMz, ip.z, ip.apexHill.accurateApexRt(specTime))), specTime)
		  				  
		  sender ! Deisotoped(batchId, uniqueIsotopes)
		  context.stop(self)
  }
  
	def deisotope(hillsByMz:Array[Hill],
			targets:Seq[Target],
			hillsMzs:Array[Double],
			hillsRts:Array[Double],
			hillsScanIndices:Array[(Int, Int)]
	):(Seq[Seq[Cluster.Edge]], Seq[IsotopePattern]) = {
		
		val targetPatterns = 
			for (t <- targets) yield {
				val monoisoHills = closeHills(hillsByMz, t, hillsMzs, hillsRts)
				if (monoisoHills.nonEmpty) {
					val patterns = getPatterns(monoisoHills, hillsByMz, t, hillsMzs, hillsScanIndices).map(ipi =>
						IsotopePattern(ipi.inds.map(hillsByMz), ipi.offset, ipi.mostAbundNbr - ipi.offset, ipi.z, ipi.averagineCorr))
					
					(t, patterns)
				} else
					(t, Nil)
			}
		
		(Nil, targetPatterns.flatMap(_._2))
	}
	
	
	def closeHills(hills:Array[Hill], t:Target, hillsMzs:Array[Double], hillsRts:Array[Double]):Seq[Int] = {
	  val (hillsStartIndx, hillsEndIndx) = DinoUtil.getMinMaxIndx(hillsMzs, t.mz, t.mzDiff)
		val inds = new ArrayBuffer[Int]
		for (i <- hillsStartIndx until hillsEndIndx) {
			if (hillsRts(i) > t.rtStart && hillsRts(i) < t.rtEnd)
				inds += i
		}
		inds
	}
	
	
	
	def getPatterns(
			seeds:Seq[Int], 
			hills:Array[Hill],
			t:Target,
			hillsMzs:Array[Double],
			hillsScanIndices:Array[(Int, Int)]
	):Seq[IsotopePatternInds] = {
		if (seeds.isEmpty) return Nil
		val patterns = 
			for {
				seed <- seeds
				isoPatInds <- extendSeed(seed, hills, t.z, hillsMzs, hillsScanIndices)
			} yield (isoPatInds)
		if (patterns.isEmpty) Nil
		else 
			patterns.filter(ip => ip.inds.head == ip.seed)
	}
	
	
	def extendSeed(seed:Int, hills:Array[Hill], z:Int, hillsMzs:Array[Double], hillsScanIndices:Array[(Int, Int)]):Option[IsotopePatternInds] = {
		
		val seedCenterMz = hillsMzs(seed)
		val seedErr2 = hills(seed).total.centerMzError * hills(seed).total.centerMzError
		val massErrorSq = params.adv.deisoSigmas * params.adv.deisoSigmas * (
					  2.5*seedErr2)
		val err2 = DinoUtil.SULPHUR_SHIFT * DinoUtil.SULPHUR_SHIFT / (z*z) + massErrorSq
		val err = Math.sqrt(err2)
		
		def extend2(direction:Int):Seq[Int] = {
			var nIso = 1
			var m = seedCenterMz + direction * nIso * DinoUtil.ISOTOPE_PATTERN_DIFF / z
			var isoMissing = false
			val isos = new ArrayBuffer[Int]
			
			while (!isoMissing) {
			  val (minIdx, maxIdx) = DinoUtil.getMinMaxIndx(hillsMzs, m, err)
			  var maxCorr = (-1, -1.0)
			  for (ii <- minIdx until maxIdx) {
	  		  if (math.min(hillsScanIndices(ii)._2, hillsScanIndices(seed)._2) > 
	  		      math.max(hillsScanIndices(ii)._1, hillsScanIndices(seed)._1)) {
	  		    val corr = params.adv.deisoCorrCalc(hills(seed), hills(ii))
	  		    if (corr > maxCorr._2) {
	  		      maxCorr = (ii, corr)
	  		    }
	  		  }
        }
				if (maxCorr._2 >= params.adv.deisoCorr) {
					isos += maxCorr._1
					nIso += 1
					m = seedCenterMz + direction * nIso * DinoUtil.ISOTOPE_PATTERN_DIFF / z
				} else 
					isoMissing = true
			}
			isos
		}
		
		/* 
		 * Whereas in global mode we also search for downMatches (instead of 
		 * just upMatches), in targeted mode this could cause the reported 
		 * precursor m/z to be lower than the one we 'target'. This lowers 
		 * the sensitivity in a matches-between-runs setting and is therefore 
		 * intentionally skipped here
		 */
	  val direction = 1
		val upMatches = extend2(direction)
		
		val result = seed +: upMatches
		val resSeedInd = 0
		val resultProfile = result.map(hills(_).apex.intensity)
		
		val minima = DinoUtil.localMinima(resultProfile, params.adv.deisoValleyFactor)
		val oneMaxResult = 
			if (minima.nonEmpty) {
				val lower = minima.filter(_ < resSeedInd).lastOption.getOrElse(0)
				val upper = minima.filter(_ > resSeedInd).headOption.getOrElse(result.length)
				result.slice(lower, upper + 1)
			} else result
		
		val cleanResult =
			if (z * seedCenterMz < 1000) {
				val apex = oneMaxResult.maxBy(hills(_).total.intensity)
				oneMaxResult.drop(oneMaxResult.indexOf(apex))
			} else oneMaxResult
		
		val cleanProfile = cleanResult.map(hills(_).total.intensity)
		val avgIsotopeDistr = Peptide.averagine((seedCenterMz - Constants.PROTON_WEIGHT)*z).getIsotopeDistribution()
		
		val alignment = params.adv.deisoAveragineStrategy(cleanProfile, avgIsotopeDistr, params)
		
		/* 
		 * In global mode we can drop some lower m/z isotopes if this causes 
		 * better agreement with the averagine isotope pattern. In a 
		 * matches-between-runs setting this again lowers the sensitivity 
		 * and is therefore intentionally skipped
		 */
		alignment.map(a => IsotopePatternInds(
				cleanResult,
				0,
				0,
				z, 
				a.corr,
				seed
			))
	}
}

case class BatchDeisotope(batchId:Int, hillsByMz:Array[Hill], batchTargets:Seq[Target], hillsMzs:Array[Double], hillsRts:Array[Double], hillsScanIndices:Array[(Int, Int)], specTime:Seq[Double])
case class Deisotoped(batchId:Int, isotopes:Seq[IsotopePattern])
case class TargetDeisotope(hills:Seq[Hill], targets:Seq[Target], specTime:Seq[Double])
case class TargetDeisotopingComplete(clusters:Seq[Seq[Cluster.Edge]], patterns:Seq[IsotopePattern])
