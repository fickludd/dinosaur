package se.lth.immun

object Estimation {
	trait EstType
	case class GaussianEst() extends EstType
	case class WeightedEst() extends EstType
	
	def center(
			minInd:Int, 
			centerInd:Int, 
			maxInd:Int, 
			xs:Seq[Double],
			ys:Seq[Double]
	)(
		implicit params:DinosaurParams	
	):(Double, Double) = {
		
		val peakIntensity =
			if (params.adv.maxIntensity)
				ys.slice(minInd, maxInd).max
			else 
				ys.slice(minInd, maxInd).sum
		
		if (minInd == maxInd) 
			return (xs(maxInd), peakIntensity)
		
		if (minInd == centerInd)
			return (est2(
						xs(centerInd), xs(centerInd + 1), 
						ys(centerInd), ys(centerInd + 1)), 
					peakIntensity)
		
		if (maxInd == centerInd) 
			return (est2(
						xs(centerInd-1), xs(centerInd), 
						ys(centerInd-1), ys(centerInd)), 
					peakIntensity)
		
		if (params.adv.massEstPoints <= 3) {
			return (
				params.adv.massEstType match {
					case GaussianEst() =>
						est3(
							xs(centerInd - 1), xs(centerInd), xs(centerInd + 1),
					        ys(centerInd - 1), ys(centerInd), ys(centerInd + 1))
					case WeightedEst() =>
						weightedMean(
							xs.slice(centerInd - 1, centerInd + 2),
					        ys.slice(centerInd - 1, centerInd + 2))
				}, peakIntensity)
		}
		
		// params.massEstPoints > 3
		var nleft = 0
		var nright = 0
		if (params.adv.massEstPoints % 2 == 1) {
			val d = params.adv.massEstPoints / 2
			nleft = math.max(centerInd - d, minInd)
			nright = math.min(centerInd + d, maxInd)
		} else {
			val d = params.adv.massEstPoints / 2 - 1
			nleft = math.max(centerInd - d, minInd)
			nright = math.min(centerInd + d, maxInd)
			if (nleft != minInd && nright != maxInd) {
				if (ys(nleft - 1) > ys(nright + 1)) 
					nleft -= 1
				else 
					nright += 1
				
			} 
			else if (nleft != minInd) nleft -= 1
			else if (nright != maxInd) nright += 1
		}
		return (estN(xs.slice(nleft, nright+1), ys.slice(nleft, nright+1), centerInd - nleft), peakIntensity)
	}
	
	def est2(mz1:Double, mz2:Double, int1:Double, int2:Double) = 
		(mz1 * int1 + mz2 * int2) / (int1 + int2)
		
	def est3(m1:Double, m2:Double, m3:Double, i1:Double, i2:Double, i3:Double) = {
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
}