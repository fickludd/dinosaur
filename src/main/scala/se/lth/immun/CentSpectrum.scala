package se.lth.immun

import Dinosaur._

class CentSpectrum(
		val index:Int,
		val minMz:Double,
		val maxMz:Double,
		val cpeaks:Seq[CentroidPeak]
) {
	
	def closestWithin(mz:Double, ppm:Double):Option[CentroidPeak] = {
		if (mz < minMz || mz > maxMz) 
			return None
		val closest = cpeaks.minBy(cp => math.abs(cp.mz - mz))
		if (2 * 1000000 * math.abs(closest.mz - mz) / (closest.mz + mz) < ppm) 
			Some(closest)
		else None
	}
	
}