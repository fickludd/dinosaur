package se.lth.immun

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.HashMap


import se.lth.immun.chem.Peptide
import se.lth.immun.chem.Constants
import se.lth.immun.chem.IsotopeDistribution

object Cluster {
	
	case class Edge(i:Int, j:Int, z:Int, corr:Double)
	
	val NO_ISOTOPE_PATTERN_INDS = IsotopePatternInds(Nil, 0, 0, 0, 0.0, 0)
	case class IsotopePatternInds(
			val inds:Seq[Int], val offset:Int, val mostAbundNbr:Int,
			val z:Int, val averagineCorr:Double, seed:Int) {
		def length = inds.length
	}
	
	case class AveragineAlignment()
	
	def apply(edges:Seq[Edge], hills:Seq[Hill])(implicit params:DinosaurParams):Cluster = {
		val indSeq = (edges.map(_.i) ++ edges.map(_.j)).distinct
		val inds = indSeq.sortBy(i => hills(i).total.centerMz)
		val l = inds.size
		val corrs = new HashMap[(Int, Int), Double]
		for (e <- edges)
			corrs += (inds.indexOf(e.j), inds.indexOf(e.i)) -> e.corr
		new Cluster(inds.map(hills), corrs)
	}
	
	
	
	def averagineFullCorr(
			profile:Seq[Double], 
			averagine:IsotopeDistribution, 
			params:DinosaurParams
	):Option[DinoUtil.Alignment] = {
		val align = DinoUtil.optimalAlign(profile, averagine.intensities)
		
		if (align.corr < params.adv.averagineCorr)
			None
		else
			Some(align)
		
		
	}
	
	
	
	def averagineOverlap(
			profile:Seq[Double], 
			averagine:IsotopeDistribution, 
			params:DinosaurParams
	):Option[DinoUtil.Alignment] = {
		val align = DinoUtil.optimalAlignOverlap(profile, averagine.intensities)
		
		if (align.corr < params.adv.averagineCorr || align.explained < params.adv.averagineExplained)
			None
		else
			Some(align)
	}
}

/*
 * Hills are sorted by mz
 */
class Cluster(
		val hills:Seq[Hill],
		val corrs:HashMap[(Int, Int), Double]
)(implicit val params:DinosaurParams) {
	
	import Cluster._
	
	val masses = hills.map(_.total.centerMz)
	val errs = hills.map(_.total.centerMzError)
	val intensities = hills.map(_.total.intensity)
	val apexIntensities = hills.map(_.apex.intensity)
	
	
	
	def deconvolve(
			inds:Seq[Int] = 0 until hills.length
	):List[IsotopePattern] = { 
		if (inds.map(hills).exists(h => params.closeToDebug(h.total.centerMz)))
			{ val k = 1 }
		getBestPattern(inds) match {
			case Some((ip, indsLeft)) =>
				IsotopePattern(ip.inds.map(hills), ip.offset, ip.mostAbundNbr - ip.offset, ip.z, ip.averagineCorr) :: deconvolve(indsLeft)
			case None =>
				Nil
		}
	}
	
	
	
	def getBestPattern(inds:Seq[Int]):Option[(IsotopePatternInds, Seq[Int])] = {
		if (inds.isEmpty) return None
		val sorted = inds.sortBy(i => hills(i).total.centerMz)
		val patterns = 
			for {
				z <- params.minCharge.value to params.maxCharge
				seed <- inds.sortBy(i => hills(i).total.intensity).reverse.take(params.adv.deisoDecompMaxSeeds)
				isopatInds <- extendSeed(seed, sorted, z)
			} yield (isopatInds)
		if (patterns.isEmpty) None
		else {
			val maxLength = patterns.map(_.length).max
			val longestPattern = patterns.filter(_.length == maxLength).maxBy(p => p.inds.map(i => hills(i).total.intensity).sum)
			Some((
					longestPattern, 
					inds.filter(i => !longestPattern.inds.contains(i))
				))
		}
	}
	
	
	
	def extendSeed(seed:Int, inds:Seq[Int], z:Int):Option[IsotopePatternInds] = {
		
		def extend2(dir:Int):Seq[Int] = {
			var ii = inds.indexOf(seed)+dir
			var nIso = 1
			var m = masses(seed) + dir * nIso * DinoUtil.ISOTOPE_PATTERN_DIFF / z
			var isoMissing = false
			val isos = new ArrayBuffer[Int]
			val alts = new ArrayBuffer[Int]
			
			def evalAlts = {
				if (alts.nonEmpty) {
					val maxCorr = alts.maxBy(a => corr(seed, a))
					alts.clear
					if (corr(seed, maxCorr) >= params.adv.deisoCorr) {
						isos += maxCorr
						nIso += 1
						m = masses(seed) + dir * nIso * DinoUtil.ISOTOPE_PATTERN_DIFF / z
					} else 
						isoMissing = true
				} else
					isoMissing = true
			}
			
			while (!isoMissing && ii < inds.length && ii >= 0) {
				val i = inds(ii)
				val massErrorSq = params.adv.deisoSigmas * params.adv.deisoSigmas * (
						errs(seed) * errs(seed) + errs(i) * errs(i))
				val err2 = DinoUtil.SULPHUR_SHIFT * DinoUtil.SULPHUR_SHIFT / (z*z) + massErrorSq
				val mDiff = masses(i) - m
				if (mDiff * mDiff <= err2) {
					alts += i
					ii += dir
				} else if (dir*mDiff > 0) {
					evalAlts
				} else
					ii += dir
			}
			evalAlts
			isos
		}
					
		val upMatches = extend2(1)
		val downMatches = extend2(-1)
		val result = downMatches.reverse ++ (seed +: upMatches)
		val resSeedInd = downMatches.length
		val resultProfile = result.map(apexIntensities)
		if (DinoUtil.isLocalMinimum(resultProfile)(resSeedInd))
			return None
		
		val minima = DinoUtil.localMinima(resultProfile, params.adv.deisoValleyFactor)
		val oneMaxResult = 
			if (minima.nonEmpty) {
				val lower = minima.filter(_ < resSeedInd).lastOption.getOrElse(0)
				val upper = minima.filter(_ > resSeedInd).headOption.getOrElse(result.length)
				result.slice(lower, upper + 1)
			} else result
					
		val cleanResult =
			if (z * masses(seed) < 1000) {
				val apex = oneMaxResult.maxBy(intensities)
				oneMaxResult.drop(oneMaxResult.indexOf(apex))
			} else oneMaxResult
		
		val cleanProfile = cleanResult.map(intensities)
		val avgIsotopeDistr = Peptide.averagine((masses(seed)-Constants.PROTON_WEIGHT)*z).getIsotopeDistribution()
		
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
	
	
	
	def corr(i1:Int, i2:Int) = {
		val min = math.min(i1, i2)
		val max = math.max(i1, i2)
		if (!corrs.contains(min, max))
			corrs += (min, max) -> 
				params.adv.deisoCorrCalc(hills(i1), hills(i2))
		corrs(min, max)
	}
}