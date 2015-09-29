package se.lth.immun

import scala.collection.mutable.ArrayBuffer

import Deisotoper._
import DinoUtil._

case class ChargePair(i1:IsotopePattern, i2:IsotopePattern)

class ChargePairer(val params:DinosaurParams) {

	def pairPatterns(patterns:Seq[IsotopePattern]):Seq[ChargePair] = {
		val sorted = patterns.sortBy(_.hills.map(_.scanIndex.last).max).toArray
		
		var minInd = 0
		val pairs = new ArrayBuffer[ChargePair]
		for (i <- 0 until sorted.length) {
			
			val pi = sorted(i)
			val minScanIndex = pi.startInd
			while (sorted(minInd).endInd < minScanIndex)
				minInd += 1
			for (j <- minInd until i) {
				val pj = sorted(j)
				
				if (within(pi.mass, pj.mass, params.adv.chargePairPPM.value)) {
					val corr = pi.totalProfile.correlate(pj.totalProfile)
					if (corr > params.adv.chargePairCorr) 
						pairs += ChargePair(pi, pj)
						
				}
			}
		}
		pairs
	}
}