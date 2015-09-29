package se.lth.immun

import scala.collection.mutable.ArrayBuffer

import se.lth.immun.chem.Peptide
import se.lth.immun.chem.Constants
import se.lth.immun.chem.IsotopeDistribution

class TargetDeisotoper(val params:DinosaurParams) extends Timeable {

	import Cluster._
	
	implicit val p = params
	
	def deisotope(
			hills:Array[Hill], 
			targets:Seq[Target], 
			specTime:Seq[Double]
	):(Seq[Seq[Cluster.Edge]], Seq[IsotopePattern]) = {
		
		println("=== IN TARGETED MODE - BEWARE!!! ===")
		val hillsByMz = hills.sortBy(_.total.centerMz)
		
		if (params.verbose)
			println("deisotoping based on targets...")
		
		val targetPatterns = 
			for (t <- targets) yield {
				if (params.verbose)
					println(t.id)
				
				val monoisoHills = closeHills(hillsByMz, t, specTime)
				if (monoisoHills.nonEmpty) {
					val patterns = getPatterns(monoisoHills, hillsByMz, t).map(ip =>
						IsotopePattern(ip.inds.map(hillsByMz), ip.offset, ip.mostAbundNbr - ip.offset, ip.z, ip.averagineCorr))
					
					(t, patterns)
				} else
					(t, Nil)
			}
		
		if (params.verbose)
			println("deisotoping complete")
		(Nil, targetPatterns.flatMap(_._2))
	}
	
	
	def closeHills(hills:Array[Hill], t:Target, specTime:Seq[Double]):Seq[Int] = {
		val inds = new ArrayBuffer[Int]
		for (i <- 0 until hills.length) {
			val h = hills(i)
			val hmz = h.total.centerMz
			val hApexRt = h.accurateApexRt(specTime)
			if (hmz > t.mz - t.mzDiff && hmz < t.mz + t.mzDiff && hApexRt > t.rtStart && hApexRt < t.rtEnd)
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
				isopatInds <- extendSeed(seed, overlapping(seed, hills), hills, t.z)
			} yield (isopatInds)
		if (patterns.isEmpty) Nil
		else 
			patterns.filter(ip => ip.inds.head == ip.seed)
	}
	
	
	def extendSeed(seed:Int, inds:Seq[Int], hills:Array[Hill], z:Int):Option[IsotopePatternInds] = {
		
		val seedTot = hills(seed).total
			
		def extend2(dir:Int):Seq[Int] = {
			var ii = inds.indexOf(seed)+dir
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
			
			while (!isoMissing && ii < inds.length && ii >= 0) {
				val iHill = hills(inds(ii))
				val iTot = iHill.total
				val mDiff = iTot.centerMz - m
				val massErrorSq = params.adv.deisoSigmas * params.adv.deisoSigmas * (
						seedErr2 + iTot.centerMzError * iTot.centerMzError)
				val err2 = DinoUtil.SULPHUR_SHIFT * DinoUtil.SULPHUR_SHIFT / (z*z) + massErrorSq
				if (mDiff * mDiff <= err2) {
					alts += inds(ii)
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
		val downMatches = extend2(-1)
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
				cleanResult.drop(math.max(0, -a.offset)), 
				math.max(0, a.offset), 
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