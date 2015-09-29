package se.lth.immun

import org.apache.commons.math3.stat.StatUtils
import scala.util.Random
import java.io.File


object DinoUtil {
	val C13C12_DIFF = 1.0033548378
	val S34S32_DIFF = 33.96786683 - 31.97207069
	val SULPHUR_SHIFT = 2 * C13C12_DIFF - S34S32_DIFF
	val ISOTOPE_PATTERN_DIFF = 1.00286864
	
	def smoothMean(xs:Seq[Double], w:Int):Seq[Double] = 
		for (i <- 0 until xs.length) yield {
			val minInd = math.max(0, i - w)
			val maxInd = math.min(xs.length, i + w + 1)
			xs.slice(minInd, maxInd).sum / (maxInd - minInd)
		}
	
	def smoothMedian(xs:Seq[Double], w:Int):Seq[Double] = 
		for (i <- 0 until xs.length) yield {
			val minInd = math.max(0, i - w)
			val maxInd = math.min(xs.length, i + w + 1)
			StatUtils.percentile(xs.slice(minInd, maxInd).toArray, 50)
		}
	
	def smoothMedian1(xs:Seq[Double]):Seq[Double] = 
		for (i <- 0 until xs.length) yield {
			if (i == 0) xs(0)
			else if (i == xs.length - 1) xs(i)
			else {
				if (xs(i-1) > xs(i)) {
					if (xs(i) > xs(i+1)) xs(i)
					else math.min(xs(i-1), xs(i+1))
				} else {
					if (xs(i-1) > xs(i+1)) xs(i-1)
					else math.min(xs(i), xs(i+1))
				}
			}
		}
	
	/**
	 * Return indices of local minima, in order of appearance
	 */
	def localMinPositions(y:Seq[Double]):Seq[Int] = {
		(2 until (y.length - 2)).filter(i =>
			isMin(y(i-2), y(i-2), y(i), y(i+1), y(i+2))
		)
	}
	
	
	def within(mz1:Double, mz2:Double, ppm:Double) = 
		2 * 1000000 * math.abs(mz1 - mz2) /	(mz1 + mz2) < ppm
	
	
		
	def isMax(m2:Double, m1:Double, x:Double, p1:Double, p2:Double):Boolean = {
		if 	(x > m1 && x > p1) 				return true
		if 	(x > m2 && x == m1 && x > p1) 	return true
		if 	(x > m1 && x == p1 && x > p2) 	return true
		false
	}
	
	def isMin(m2:Double, m1:Double, x:Double, p1:Double, p2:Double):Boolean = {
		if (x < m1 && x < p1) 						return true
		if (x == m1 && x < m2 && x < p1) 			return true
		if (x < m1 && x == p1 && x < p2) 			return true
		if (x < m2 && x == m1 && x == p1 && x < p2) return true
		false
	}
	
	def isLocalMinimum(xs:Seq[Double])(i:Int):Boolean = {
		if (i == 0 || i == xs.length - 1) false
		else
			xs(i-1) > xs(i) && xs(i+1) > xs(i)
	}
	
	def localMinima(xs:Seq[Double], valleyFactor:Double):Seq[Int] = {
		val mins = (1 until xs.length-1).filter(DinoUtil.isLocalMinimum(xs))
		
		mins.filter(min => {
			val leftMax = xs.take(min).max
			val rightMax = xs.drop(min+1).max
			val smallMax = math.min(leftMax, rightMax)
			smallMax / xs(min) >= valleyFactor
		})
	}
	
	def ensureDir(dir:File) = {
		if (!dir.exists)
			dir.mkdir
	}
	

		
	def cosineCorrelation(x:Seq[Double], y:Seq[Double]):Double = {
		if (x.length < 3) return 0
		
		var xx = 0.0
		var yy = 0.0
		var xy = 0.0
		var i = 0
		while (i < x.length) {
			val wx = x(i)
			val wy = y(i)
			xx += wx * wx
			yy += wy * wy
			xy += wx * wy
			i += 1
		}
		val denom = xx * yy
		if (denom > 0.0) 
			xy / math.sqrt(denom)
		else 0.0
	}
	
	def cosineCorrelation(p:Seq[(Double,Double)]):Double = {
		if (p.length < 3) return 0
		
		var xx = 0.0
		var yy = 0.0
		var xy = 0.0
		var i = 0
		while (i < p.length) {
			val (wx, wy) = p(i)
			xx += wx * wx
			yy += wy * wy
			xy += wx * wy
			i += 1
		}
		val denom = xx * yy
		if (denom > 0.0) 
			xy / math.sqrt(denom)
		else 0.0
	}
	
	
	
	def cosineCorrelationSparse(ixs:Seq[Int], xs:Seq[Double], iys:Seq[Int], ys:Seq[Double])(implicit params:DinosaurParams):Double = {
		var ix = 0
		var iy = 0
		var xx = 0.0
		var yy = 0.0
		var xy = 0.0
		var overlap = 0
		while (ix < ixs.length && iy < iys.length) {
			if (ixs(ix) == iys(iy)) {
				val wx = xs(ix)
				val wy = ys(iy)
				xx += wx * wx
				xy += wx * wy
				yy += wy * wy
				ix += 1
				iy += 1
				overlap += 1
			} else {
				if (ixs(ix) < iys(iy)) 
					ix += 1
				else 
					iy += 1
			}
		}
		
		val denom = xx * yy
		if (denom > 0.0 && overlap >= params.adv.deisoOverlap) 
			xy / math.sqrt(denom)
		else 0.0
	}
	
	
	
	def cosineCorrelationFull(ixs:Seq[Int], xs:Seq[Double], iys:Seq[Int], ys:Seq[Double]):Double = {
		var ix = 0
		var iy = 0
		var xx = 0.0
		var yy = 0.0
		var xy = 0.0
		var overlap = 0
		while (ix < ixs.length && iy < iys.length) {
			val wx = xs(ix)
			val wy = ys(iy)
			val scanX = ixs(ix)
			val scanY = iys(iy)
			if (scanX == scanY) {
				xy += wx * wy
				overlap += 1
			} 
			if (scanX <= scanY) {
				xx += wx * wx
				ix += 1
			} 
			if (scanX >= scanY) {
				yy += wy * wy
				iy += 1
			}
		}
		
		while (ix < xs.length) {
			val wx = xs(ix)
			xx += wx * wx
			ix += 1
		}
		while (iy < ys.length) {
			val wy = ys(iy)
			yy += wy * wy
			iy += 1
		}
		
		val denom = xx * yy
		if (denom > 0.0) 
			xy / math.sqrt(denom)
		else 0.0
	}
	
	
	
	
	
	
	case class Alignment(offset:Int, corr:Double, explained:Double)
	
	/**
	 * Return pair with offset and correlation to the reference
	 * The offset indicates how many steps xs should be shifted 
	 * to the right for the highest correlation between the two
	 */
	
	
	def optimalAlign(xs:Seq[Double], ref:Seq[Double]):Alignment = {
		val yy = ref.map(y => y*y).sum
		val offCorrs = 
			for (off <- -1 until xs.length - 1) yield {
				var xx = 0.0
				var xy = 0.0
				var explainedMass = 0.0
				for (i <- 0 until ref.length) {
					val ix = i-off
					if (ix >= 0 && ix < xs.length) {
						val y = ref(i)
						val x = xs(ix)
						xx += x * x
						xy += x * y
						explainedMass += y
					}
				}
				val denom = xx * yy
				if (denom > 0.0) 
					Alignment(off, xy / math.sqrt(denom), explainedMass)
				else Alignment(off, 0.0, explainedMass)
			}
		offCorrs.maxBy(_.corr)
	}
	
	
	def optimalAlignOverlap(xs:Seq[Double], ref:Seq[Double]):Alignment = {
		if (xs.length <= 1)
			return Alignment(0, 0.0, 0.0)
		val offCorrs = 
			for (off <- 2 - xs.length to xs.length - 2) yield {
				var xx = 0.0
				var xy = 0.0
				var yy = 0.0
				var explainedMass = 0.0
				for (i <- 0 until ref.length) {
					val ix = i-off
					if (ix >= 0 && ix < xs.length) {
						val y = ref(i)
						val x = xs(ix)
						xx += x * x
						xy += x * y
						yy += y * y
						explainedMass += y
					}
				}
				val denom = xx * yy
				if (denom > 0.0) 
					Alignment(off, xy / math.sqrt(denom), explainedMass)
				else Alignment(off, 0.0, explainedMass)
			}
		offCorrs.maxBy(t => t.corr * t.explained)
	}
	
	
	def alignApex(x1:Seq[Double], x2:Seq[Double]):(Seq[Double], Seq[Double]) = {
		val im1 = (0 until x1.length).maxBy(x1)
		val im2 = (0 until x2.length).maxBy(x2)
		
		def pad(s:Seq[Double], offset:Int) =
			s.padTo(math.max(x1.length - im1, x2.length - im2)+offset, 0.0)
		
		if (im1 == im2) (x1, x2)
		else if (im1 < im2) 	(pad(new Array[Double](im2 - im1) ++ x1, im2), pad(x2, im2))
		else 					(pad(x1, im1), pad(new Array[Double](im1 - im2) ++ x2, im1))
	}
	
	
	def interpolateLinear(x:Double, x0:Double, x1:Double, y0:Double, y1:Double):Double = {
		val k = (x-x0) / (x1-x0)
		k*y1 + (1-k)*y0
	}
}