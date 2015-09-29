package se.lth.immun

import org.apache.commons.math3.analysis.differentiation.MultivariateDifferentiableVectorFunction
import org.apache.commons.math3.optimization.general.LevenbergMarquardtOptimizer
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure
import org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction
import org.apache.commons.math3.analysis.MultivariateMatrixFunction


case class MzIntensity(mz:Double, intensity:Double)

object ChargePairs {
	val m0 = 445.120025
	val mConst = 100/math.sqrt(m0)
	val intConst = math.log(1e7)
}

class ChargePairsDiff(pairs:Seq[ChargePair]) extends DifferentiableMultivariateVectorFunction {
	
	
	
	def value(variables:Array[Double]):Array[Double] = {
		(for (ChargePair(ip1, ip2) <- pairs) yield {
			f(ip1.mass, ip1.intensity, ip2.mass, ip2.intensity)(variables)
		}).toArray
	}
	
	
	
	def target:Array[Double] = pairs.map(_ => 0.0).toArray
	
	
	def weights = 
		pairs.map(_ match {
			case ChargePair(ip1, ip2) => 
				ip1.massError * ip1.massError + ip2.massError * ip2.massError
		}).toArray
	
	
	def jacobian():MultivariateMatrixFunction = 
		new MultivariateMatrixFunction {
			def value(point:Array[Double]) = 
                jacobian(point)
		}
	
	def jacobian(variables:Array[Double]):Array[Array[Double]] = 
		(for (ChargePair(ip1, ip2) <- pairs) yield {
			dfs(ip1.mass, ip1.intensity, ip2.mass, ip2.intensity, variables)
		}).toArray
	
		
		
	import ChargePairs._
	def fPart(mz:Double, intensity:Double, pq:Seq[Double]):Double = {
		val u = 100/math.sqrt(mz) - mConst
		val v = math.log(intensity) - intConst
		val P = u*(pq(0) + u*(pq(1) + u*(pq(2) + u*(pq(3) + u*pq(4)))))
		val Q = pq(5)*v
		mz * (1 + 1e-6*(P + Q))
	}
	
	
	
	def f(mz1:Double, intensity1:Double, mz2:Double, intensity2:Double)(pq:Seq[Double]):Double = 
		fPart(mz1, intensity1, pq) - fPart(mz2, intensity2, pq)
	
	
	
	val epsilon = 1e-7
	def dfs(mz1:Double, int1:Double, mz2:Double, int2:Double, pq:Array[Double]):Array[Double] = {
		val nd = numDiff(mz1, int1, mz2, int2, pq) _
		Array(nd(0, epsilon), nd(1, epsilon), nd(2, epsilon), nd(3, epsilon), nd(4, epsilon), nd(5, epsilon))
	}
	
	
	
	def dPQ(pq:Seq[Double], i:Int, diff:Double):Seq[Double] =
			pq.updated(i, pq(i)+diff)
	
	
	
	def numDiff(mz1:Double, int1:Double, mz2:Double, int2:Double, pq:Array[Double])(i:Int, epsilon:Double):Double = {
		val g = f(mz1, int1, mz2, int2) _
		(g(dPQ(pq, i, epsilon)) - g(dPQ(pq, i, -epsilon))) / (2*epsilon)
	}
	
}


class ChargePairsOne(pairs:Seq[ChargePair]) extends DifferentiableMultivariateVectorFunction {
	
	val xs = pairs.map(_.i1) ++ pairs.map(_.i2)
	val ys = pairs.map(_.i2.mass) ++ pairs.map(_.i1.mass)
	
	def value(variables:Array[Double]):Array[Double] = {
		(for (ip <- xs) yield {
			f(ip.mass, ip.intensity)(variables)
		}).toArray
	}
	
	
	
	def target:Array[Double] = ys.toArray
	
	
	def weights:Array[Double] = {
		val w = pairs.map(_ match {
			case ChargePair(ip1, ip2) => 
				ip1.massError * ip1.massError + ip2.massError * ip2.massError
		}).toArray
		(w ++ w).toArray
	}
		
	
	def jacobian():MultivariateMatrixFunction = 
		new MultivariateMatrixFunction {
			def value(point:Array[Double]) = 
                jacobian(point)
		}
	
	def jacobian(variables:Array[Double]):Array[Array[Double]] = 
		(for (ip <- xs) yield {
			dfs(ip.mass, ip.intensity, variables)
		}).toArray
	
		
		
	import ChargePairs._
	def f(mz:Double, intensity:Double)(pq:Seq[Double]):Double = {
		val u = 100/math.sqrt(mz) - mConst
		val v = math.log(intensity) - intConst
		val P = u*(pq(0) + u*(pq(1) + u*(pq(2) + u*(pq(3) + u*pq(4)))))
		val Q = pq(5)*v
		mz * (1 + 1e-6*(P + Q))
	}
	
	
	
	val epsilon = 1e-7
	def dfs(mz:Double, int:Double, pq:Array[Double]):Array[Double] = {
		val nd = numDiff(mz, int, pq) _
		Array(nd(0, epsilon), nd(1, epsilon), nd(2, epsilon), nd(3, epsilon), nd(4, epsilon), nd(5, epsilon))
	}
	
	
	
	def dPQ(pq:Seq[Double], i:Int, diff:Double):Seq[Double] =
			pq.updated(i, pq(i)+diff)
	
	
	
	def numDiff(mz:Double, int:Double, pq:Array[Double])(i:Int, epsilon:Double):Double = {
		val g = f(mz, int) _
		(g(dPQ(pq, i, epsilon)) - g(dPQ(pq, i, -epsilon))) / (2*epsilon)
	}
	
}


class NonLinMassCalibration(val params:DinosaurParams) {
	
	import ChargePairs._
	
	def findRecalibrationFunction(pairs:Seq[ChargePair]):MzIntensity => Double = {
		val problem = new ChargePairsOne(pairs)
		val optimizer = new LevenbergMarquardtOptimizer
		val initialSolution = Array(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
		val optimum = optimizer.optimize(100, problem, problem.target, problem.weights, initialSolution)
		val pq = optimum.getPoint
		_ match {
			case MzIntensity(mz, intensity) =>
				val u = 100/math.sqrt(mz) - mConst
				val v = math.log(intensity) - intConst
				val P = u*(pq(0) + u*(pq(1) + u*(pq(2) + u*(pq(3) + u*pq(4)))))
				val Q = pq(5)*v
				mz * (1 + 1e-6*(P + Q))
		}
	}
}