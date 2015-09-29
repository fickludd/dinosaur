package se.lth.immun

import scala.util.Random

object Bootstrap {
	val bufferLen = 999
	val maxVectorLen = 99
	val buffer = {
		val a = new Array[Array[Array[Int]]](bufferLen)
		for (i <- 0 until bufferLen) {
			a(i) = new Array[Array[Int]](maxVectorLen)
			for (j <- 0 until maxVectorLen) {
				a(i)(j) = new Array[Int](j)
				for (k <- 0 until j)
					a(i)(j)(k) = Random.nextInt(j)
			}
		}
		a
	}
	var count = 0
	
	
	def set(k:Int, nBoot:Int):Seq[Int] = set(k, k, nBoot)
	
	def set(k:Int, n:Int, nBoot:Int):Seq[Int] = 
		if (k<=n && n < maxVectorLen && nBoot < bufferLen) {
			count = (count + 1) % bufferLen
			if (k == n)
				buffer(count)(n)
			else
				buffer(count)(n).take(k)
		} else 
			(0 until k).map(_ => Random.nextInt(n))
	
	val table2total = 4
	val table2 = Array(
		Array(2,0,1),
		Array(1,1,2),
		Array(0,2,1)
	)
	
	val table3total = 27
	val table3 = Array(
		Array(3,0,0,1),
		Array(2,1,0,3),
		Array(2,0,1,3),
		
		Array(0,3,0,1),
		Array(1,2,0,3),
		Array(0,2,1,3),
		
		Array(0,0,3,1),
		Array(1,0,2,3),
		Array(0,1,2,3),
		
		Array(1,1,1,6)
	)
	
	
	def weightedAverage(xs:Seq[Double], ws:Seq[Double], nBoot:Int, maxBootSize:Int):(Double, Double) = {
		val l = xs.length
		val wXs = new Array[Double](l)
		for (i <- 0 until l)
			wXs(i) = ws(i) * xs(i)
		
		if (l > 3)
			bootstrapWeightedAverage(wXs, ws, nBoot, maxBootSize)
		else if (l == 2)
			weightedAverageFromTable(Bootstrap.table2, wXs, ws)
		else if (l == 3) 
			weightedAverageFromTable(Bootstrap.table3, wXs, ws)
		else
			throw new Exception("This should never happen!")
	}
	
	def weightedAverageFromTable(
			tab:Array[Array[Int]], 
			wXs:Seq[Double], 
			ws:Seq[Double]
	):(Double, Double) = {
		val tl = tab.length
		val ys = new Array[Double](tl)
		val ks = new Array[Int](tl)
		var sumYk = 0.0
		var sumK = 0
		var i = 0
		while (i < tl) {
			val row = tab(i)
			val k = row.last
			val y = weightedAverageTable(wXs, ws, row)
			ys(i) = y
			ks(i) = k
			sumYk += y * k
			sumK += k
			i += 1
		}
		val avgY = sumYk / sumK
		var err = 0.0
		i = 0
		while (i < tl) {
			val residual = avgY - ys(i)
			err += residual * residual * ks(i)
			i += 1
		}
		(
			avgY,
			math.sqrt(err / (sumK-1))
		)
	}
	
	def weightedAverageTable(wXs:Seq[Double], ws:Seq[Double], row:Array[Int]):Double = {
		var m = 0.0
		var norm = 0.0
		var j = 0
		while (j < row.length - 1) {
			val k = row(j)
			m += wXs(j)*k
			norm += ws(j)*k
			j += 1
		}
		m / norm
	}
	
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
			val y = weightedAverage(wXs, ws, Bootstrap.set(math.min(maxBootSize, l), l, nBoot))
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