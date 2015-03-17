package se.lth.immun

import akka.actor._
import se.lth.immun.mzml.ghost.GhostSpectrum

import scala.collection.mutable.ArrayBuffer


object CentroidWorker {
	case class CentroidSpectrum(gs:GhostSpectrum, customer:ActorRef, writeReport:Boolean)
	case class CentroidingComplete(specIndex:Int, centroidedPeaks:Seq[CentroidPeak], t:Long)
	
	
	case class CentroidPeak(mz:Double, int:Double, minMz:Double, maxMz:Double)
	
	trait MassEstType
	case class GaussianMassEst() extends MassEstType
	case class WeightedMassEst() extends MassEstType
}

class CentroidWorker(val params:DinosaurParams) extends Actor {

	import CentroidWorker._
	
	def receive = {
		case CentroidSpectrum(gs, customer, writeReport) =>
			centroidSpectrum(gs, customer, writeReport)
	}
	
	def centroidSpectrum(gs:GhostSpectrum, customer:ActorRef, writeReport:Boolean) = {
		
		val t0 = System.currentTimeMillis
		
		val (mzs, ints) = filter(gs)
			
		def isMax(x:Double, m1:Double, p1:Double, m2:Double, p2:Double) = {
			if 			(x > m1 && x > p1) 				true
			else if 	(x > m2 && x == m1 && x > p1) 	true
			else if 	(x > m1 && x == p1 && x > p2) 	true
			else 										false
		}

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
		
		for (i <- 2 until (mzs.length-2)) {
			val m2 	= ints(i - 2)
			val m1 	= ints(i - 1)
			val x 	= ints(i)
			val p1 	= ints(i + 1)
			val p2 	= ints(i + 2)
			if (x >= params.intensityThreshold) {
				if (isMax(x, m1, p1, m2, p2)) {
					var minInd = minPeakIndex(i)
					var maxInd = maxPeakIndex(i)
					if (maxInd - minInd > 2) {
						
						/*
						 * This part is wierd. Try to exclude later?
						 */
						if (maxInd > i && minInd < i) {
							maxInd -= 1
							minInd += 1
						}
						else if (maxInd > i) 
							maxInd = i + 1
						else if (minInd < i)
							minInd = i - 1
						
						val (mz, int) = massCenter(minInd, i, maxInd, mzs, ints)
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
		
		if (writeReport) {
			import CentroidReportActor._
	
			val reporter = context.actorOf(Props(new CentroidReportActor(params)))
			reporter ! WriteReport(gs.spectrum.index, mzs, ints, res)
		}
		
		customer ! CentroidingComplete(gs.spectrum.index, res, System.currentTimeMillis - t0)
	}
	
	def isNumber(d:Double) = java.lang.Double.isNaN(d) == false
	
	
	// here low quantile peaks might be filtered (see params)
	def filter(gs:GhostSpectrum) = 
		(gs.mzs, gs.intensities)
	
	
	def massCenter(
			minInd:Int, 
			centerInd:Int, 
			maxInd:Int, 
			mzs:Seq[Double],
			ints:Seq[Double]
	):(Double, Double) = {
		
		var peakIntensity = 0.0
		for (j <- minInd to maxInd)
			if (params.maxIntensity) 
				if (ints(j) > peakIntensity)
					peakIntensity = ints(j)
			else 
				peakIntensity += ints(j)
		
		if (minInd == maxInd) 
			return (mzs(maxInd), peakIntensity)
		
		if (minInd == centerInd)
			return (massEst2(
						mzs(centerInd), mzs(centerInd + 1), 
						ints(centerInd), ints(centerInd + 1)), 
					peakIntensity)
		
		if (maxInd == centerInd) 
			return (massEst2(
						mzs(centerInd-1), mzs(centerInd), 
						ints(centerInd-1), ints(centerInd)), 
					peakIntensity)
		
		if (params.massEstPoints <= 3) {
			return (
				params.massEstType match {
					case GaussianMassEst() =>
						massEst3(
							mzs(centerInd - 1), mzs(centerInd), mzs(centerInd + 1),
					        ints(centerInd - 1), ints(centerInd), ints(centerInd + 1))
					case WeightedMassEst() =>
						weightedMean(
							mzs.slice(centerInd - 1, centerInd + 1),
					        ints.slice(centerInd - 1, centerInd + 1))
				}, peakIntensity)
		}
		
		// params.massEstPoints > 3
		var nleft = 0
		var nright = 0
		if (params.massEstPoints % 2 == 1) {
			val d = params.massEstPoints / 2
			nleft = math.max(centerInd - d, minInd)
			nright = math.min(centerInd + d, maxInd)
		} else {
			val d = params.massEstPoints / 2 - 1
			nleft = math.max(centerInd - d, minInd)
			nright = math.min(centerInd + d, maxInd)
			if (nleft != minInd && nright != maxInd) {
				if (ints(nleft - 1) > ints(nright + 1)) 
					nleft -= 1
				else 
					nright += 1
				
			} 
			else if (nleft != minInd) nleft -= 1
			else if (nright != maxInd) nright += 1
		}
		return (estN(mzs.slice(nleft, nright+1), ints.slice(nleft, nright+1), centerInd - nleft), peakIntensity)
	}
	
	def massEst2(mz1:Double, mz2:Double, int1:Double, int2:Double) = 
		(mz1 * int1 + mz2 * int2) / (int1 + int2)
		
	def massEst3(m1:Double, m2:Double, m3:Double, i1:Double, i2:Double, i3:Double) = {
			val l1 = math.log(i1)
			val l2 = math.log(i2)
			val l3 = math.log(i3)
			0.5 * ((l1 - l2) * (m3 * m3 - m1 * m1) - (l1 - l3) * (m2 * m2 - m1 * m1)) /
			       ((l1 - l2) * (m3 - m1) - (l1 - l3) * (m2 - m1))
		}
		
	
	def weightedMean(x:Seq[Double], y:Seq[Double]):Double = {
		var m = 0.0
		var w = 0.0
		var i = 0
		while( i < x.length) {
			m += x(i) * y(i)
			w += y(i)
			i += 1
		}
		m / w
	}
	
	def estN(x:Seq[Double], y:Seq[Double], ic:Int) = {
		val dm = x(ic)
		val xs = x.map(_ - dm)
		val ys = y.map(math.log)
		throw new Exception("not yet implemented")
	}

	/*
	def Qwert(x:Double, a:Array[Double]) = {
		a(0) = 1
		a(1) = x
		a(2) = x * x
	}
	*/
}