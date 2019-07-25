package se.lth.immun

import java.util.SplittableRandom

class SplittableBootstrap(
		val rng:SplittableRandom,
		val bufferOption:Option[Array[Array[Array[Int]]]] = None
) {
	val bufferLen = 999
	val maxVectorLen = 99
	val buffer = bufferOption getOrElse {
		val a = new Array[Array[Array[Int]]](bufferLen)
		for (i <- 0 until bufferLen) {
			a(i) = new Array[Array[Int]](maxVectorLen)
			for (j <- 0 until maxVectorLen) {
				a(i)(j) = new Array[Int](j)
				for (k <- 0 until j)
					a(i)(j)(k) = rng.nextInt(j)
			}
		}
		a
	}
	var count = 0
	
	
	def split():SplittableBootstrap = {
		val rng_splitted = rng.split()
		new SplittableBootstrap(rng_splitted, Some(buffer))
	}
	
	def set(k:Int, n:Int, nBoot:Int):Seq[Int] = 
		if (k<=n && n < maxVectorLen && nBoot < bufferLen) {
			count = (count + 1) % bufferLen
			if (k == n)
				buffer(count)(n)
			else
				buffer(count)(n).take(k)
		} else 
			(0 until k).map(_ => rng.nextInt(n))
	
	
	def bootstrapWeightedAverage(
			wXs:Seq[Double], 
			ws:Seq[Double],
			nBoot:Int,
			maxBootSize:Int
	):(Double, Double) = {
		val l = wXs.length
		val ys = new Array[Double](nBoot)
		var i = 0
		var sumY = 0.0
		while (i < nBoot) {
			val y = weightedAverage(wXs, ws, set(math.min(maxBootSize, l), l, nBoot))
			ys(i) = y
			sumY += y
			i += 1
		}
		
		val avgY = sumY / nBoot
		var err = 0.0
		i = 0
		while (i < nBoot) {
			val residual = avgY - ys(i)
			err += residual * residual
			i += 1
		}
		(avgY, math.sqrt(err / (nBoot-1)))
	}
	
	def weightedAverage(wXs:Seq[Double], ws:Seq[Double], indices:Seq[Int]):Double = {
		var ii = 0
		var m = 0.0
		var norm = 0.0
		while (ii < indices.length) {
			val i = indices(ii)
			m += wXs(i)
			norm += ws(i)
			ii += 1
		}
		m / norm
	}
}
