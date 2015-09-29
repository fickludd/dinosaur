package se.lth.immun



import se.lth.immun.mzml.ghost.GhostSpectrum
import scala.collection.mutable.ArrayBuffer



object Centroider {
	
	val MS_CENTROID_SPECTRUM = "MS:1000127"
	val MS_PROFILE_SPECTRUM = "MS:1000128"
}



class Centroider(val params:DinosaurParams, streamer:ReportStreamer) {
	
	import Centroider._
	import Dinosaur._
	import DinoUtil._
	
	def centroidSpectrum(gs:GhostSpectrum, writeReport:Boolean) = {
		
		val t0 = System.currentTimeMillis
		
		val preCentroid = gs.spectrum.cvParams.exists(_.accession == MS_CENTROID_SPECTRUM)
		val preProfile = gs.spectrum.cvParams.exists(_.accession == MS_PROFILE_SPECTRUM)
		
		if (preCentroid && !preProfile) {
			(new CentSpectrum(
					gs.spectrum.index, 
					gs.mzs.min, 
					gs.mzs.max, 
					gs.mzs.zip(gs.intensities).map(t => 
						new CentroidPeak(t._1, t._2, t._1, t._1)
					)
				), 
				System.currentTimeMillis - t0
			)
		} else if (preProfile && !preCentroid) {
			val (mzs, ints) = filter(gs)
	
			def minPeakIndex(ind:Int) = {
				var i = ind
				while (i > 0 && ints(i-1) != 0 && ints(i-1) < ints(i)) 
					i -= 1
				i
			}
			def maxPeakIndex(ind:Int) = {
				var i = ind
				while (i < ints.length-1 && ints(i+1) != 0 && ints(i+1) < ints(i)) 
					i += 1
				i
			}
			
			val res = new ArrayBuffer[CentroidPeak]()
			
			if (params.atDebugTime(gs.scanStartTime)) 
				{ val k = 1 }
			
			for (i <- 2 until (mzs.length-2)) {
				val m2 	= ints(i - 2)
				val m1 	= ints(i - 1)
				val x 	= ints(i)
				val p1 	= ints(i + 1)
				val p2 	= ints(i + 2)
				if (x >= params.adv.intensityThreshold) {
					if (isMax(m2, m1, x, p1, p2)) {
						var minInd = minPeakIndex(i)
						var maxInd = maxPeakIndex(i)
						if (maxInd - minInd > 2) {
							
							/*
							 * This part is wierd. Try to exclude later?
							 * 
							 * TRYING...
							if (maxInd > i && minInd < i) {
								maxInd -= 1
								minInd += 1
							}
							else if (maxInd > i) 
								maxInd = i + 1
							else if (minInd < i)
								minInd = i - 1
							 */
							
							val (mz, int) = Estimation.center(minInd, i, maxInd, mzs, ints)(params)
							if (isNumber(mz) && isNumber(int)) {
								val minMz = 
									if (minInd > 0)
										0.5 * (mzs(minInd) + mzs(minInd - 1))
									else 
										1.5 * mzs(0) - 0.5 * mzs(1)
								
								val maxMz = 
									if (maxInd < mzs.length - 1) 
										0.5 * (mzs(maxInd) + mzs(maxInd + 1))
									else
										1.5 * mzs(maxInd) - 0.5 * mzs(maxInd - 1)
								
								res += CentroidPeak(mz, int, minMz, maxMz)
							}
						}
					}
				}
			}
			
			if (writeReport)
				CentroidReport(gs.spectrum.index, mzs, ints, res, streamer)
				
			(new CentSpectrum(gs.spectrum.index, mzs.min, mzs.max, res), System.currentTimeMillis - t0)
		} else 
			throw new Exception("Cannot tell if spectrum '%d' is centroid of profile".format(gs.spectrum.index))
	}
	
	def isNumber(d:Double) = java.lang.Double.isNaN(d) == false
	
	
	// here low quantile peaks might be filtered (see params)
	def filter(gs:GhostSpectrum) = 
		(gs.mzs, gs.intensities)
	
	

	/*
	def Qwert(x:Double, a:Array[Double]) = {
		a(0) = 1
		a(1) = x
		a(2) = x * x
	}
	*/
}