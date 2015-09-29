package se.lth.immun

import scala.collection.mutable.ArrayBuffer


case class HillSummary(
		val minMz:Double, 
		val maxMz:Double, 
		val intensity:Double, 
		val centerMz:Double, 
		val centerMzError:Double
	)

case class ScanIntensity(val scanIndex:Int, val intensity:Double)
	
case class HillProfile(i0:Int, intensity:Seq[Double]) {
	def correlate(hp:HillProfile) = {
		val i1 = List[Double]().padTo(math.max(0, i0 - hp.i0), 0.0) ++ intensity
		val i2 = List[Double]().padTo(math.max(0, hp.i0 - i0), 0.0) ++ hp.intensity
		DinoUtil.cosineCorrelation(
				i1.padTo(math.max(i1.length, i2.length), 0.0), 
				i2.padTo(math.max(i1.length, i2.length), 0.0))
	}
}

class HillHillStats(h1:Hill, h2:Hill, z:Int, params:DinosaurParams) {
	import DinoUtil._
	
	val mErrPart = h1.total.centerMzError * h1.total.centerMzError + h2.total.centerMzError * h2.total.centerMzError
	val nStds = 5
	val mError = nStds * math.sqrt(mErrPart)
	val totError = math.sqrt(SULPHUR_SHIFT * SULPHUR_SHIFT / (z*z) + mError*mError)
	val mDiff = math.abs(h1.total.centerMz - h2.total.centerMz)
	val expMDiff = ISOTOPE_PATTERN_DIFF / z
	val resMDiff = math.abs(mDiff - expMDiff)
	val corr = params.adv.deisoCorrCalc(h1, h2)
}
	
object Hill {
	val NOSummary = HillSummary(0,0,0,0,0)

	def apply(cp:CentroidPeak, ms1Index:Int)(implicit params:DinosaurParams):Hill = {
		val h = new Hill
		h.centerMz += cp.mz
		h.updateMzGuess(cp.mz)
		h.rawIntensity += cp.int
		//h.minMz += cp.minMz
		//h.maxMz += cp.maxMz
		h.scanIndex += ms1Index
		h.lastIndex = ms1Index
		h.maxMzWidth = math.max(h.maxMzWidth, cp.maxMz - cp.minMz)
		h
	}
}



class Hill(	
	val centerMz:ArrayBuffer[Double] = new ArrayBuffer[Double],
	//val minMz:ArrayBuffer[Double] = new ArrayBuffer[Double],
	//val maxMz:ArrayBuffer[Double] = new ArrayBuffer[Double],
	val rawIntensity:ArrayBuffer[Double] = new ArrayBuffer[Double],
	val scanIndex:ArrayBuffer[Int] = new ArrayBuffer[Int]
)(implicit params:DinosaurParams) {
	var smoothIntensity:Seq[Double] = Nil
	var lastIndex = 0
	var maxMzWidth = 0.0
	
	var leftSplit = false
	var rightSplit = false
	
	var total:HillSummary = Hill.NOSummary
	
	def length = centerMz.length
	
	override def toString = "Hill(n=%d,mz=%.2f,%d-%d)".format(length, if (total.centerMz > 0) total.centerMz else centerMz.head, scanIndex.head, scanIndex.last)
	
	var mzGuess:Double = 0.0
	var mzGuessN:Int = 0
	def updateMzGuess(mz:Double) = {
		mzGuess = (mzGuess * mzGuessN + mz) / (mzGuessN + 1)
		if (mzGuessN < params.adv.hillMzGuessLength - 1)
			mzGuessN += 1
	}
	
	
	
	def apex:ScanIntensity = {
		val (scan, int) = scanIndex.zip(smoothIntensity).maxBy(_._2)
		val i = scanIndex.indexOf(scan)
		val rawSlice = rawIntensity.slice(
				math.max(0, i-1),
				math.min(scanIndex.length-1, i+1))
		ScanIntensity(scan, rawSlice.sum / rawSlice.length)
	}
	
	
	
	def accurateApexRt(specTime:Seq[Double]):Double = {
		val (scan, int) = scanIndex.zip(smoothIntensity).maxBy(_._2)
		val i = scanIndex.indexOf(scan)
		val (rt, _) = Estimation.center(0, i, length-1, scanIndex.map(specTime), smoothIntensity)
		rt
	}
	
	
	def startRt(specTime:Seq[Double]):Double = {
		val rt0 = specTime(math.max(0, scanIndex.head-1))
		val rt1 = specTime(scanIndex.head)
		(rt0 + rt1) / 2
	}
	
	
	def endRt(specTime:Seq[Double]):Double = {
		val rt0 = specTime(scanIndex.last)
		val rt1 = specTime(math.min(specTime.length-1, scanIndex.last+1))
		(rt0 + rt1) / 2
	}
	
	
	def fwhm(specTime:Seq[Double]):Double = {
		val (scan, max) = scanIndex.zip(smoothIntensity).maxBy(_._2)
		val i = scanIndex.indexOf(scan)
		var i0 = i
		while (i0 > 0 && smoothIntensity(i0) > max / 2) {
			i0 -= 1
		}
		var i1 = i
		while (i1 < length-1 && smoothIntensity(i1) > max / 2) {
			i1 += 1
		}
		val rt0 = DinoUtil.interpolateLinear(
				max/2, smoothIntensity(i0), smoothIntensity(i0+1), 
				specTime(scanIndex(i0)), specTime(scanIndex(i0+1)))
		val rt1 = DinoUtil.interpolateLinear(
				max/2, smoothIntensity(i1-1), smoothIntensity(i1), 
				specTime(scanIndex(i1-1)), specTime(scanIndex(i1)))
		rt1 - rt0
	}
	
	
	
	def push(cp:CentroidPeak, ind:Int):Hill = {
		centerMz += cp.mz
		//minMz += cp.minMz
		//maxMz += cp.maxMz
		rawIntensity += cp.int
		scanIndex += ind
		lastIndex = ind
		this
	}
	
	
	
	def clear = {
		centerMz.clear
		//minMz.clear
		//maxMz.clear
		rawIntensity.clear
		scanIndex.clear
		lastIndex = 0
	}
	
	
	
	def append(h:Hill):Hill = {
		centerMz ++= h.centerMz
		//minMz ++= h.minMz
		//maxMz ++= h.maxMz
		rawIntensity ++= h.rawIntensity
		scanIndex ++= h.scanIndex
		lastIndex = h.lastIndex
		this
	}
	
	
	
	def corr(f:(Seq[Int],Seq[Double],Seq[Int],Seq[Double]) => Double)(h:Hill):Double = 
		f(scanIndex, smoothIntensity, h.scanIndex, h.smoothIntensity)
	
	
	
	def calcAverages(implicit params:DinosaurParams):Hill = {
		val totalIntensity = rawIntensity.sum
			
		val (tMz, tMzError) = calcTotalMz
		val totalMinMz = tMz - maxMzWidth/2 //minMz.min
		val totalMaxMz = tMz + maxMzWidth/2 //maxMz.max
		total = HillSummary(totalMinMz, totalMaxMz, totalIntensity, tMz, tMzError)
		this
	}
	
	
	
	def calcTotalMz(implicit params:DinosaurParams):(Double, Double) = {
		val weightedMz = new Array[Double](length)
		for (i <- 0 until length)
			weightedMz(i) = rawIntensity(i) * centerMz(i)
		
		Bootstrap.bootstrapWeightedAverage(weightedMz, rawIntensity, params.adv.hillNBoots, params.adv.maxBootSize)
	}
	
		
	
	def smooth = {
		//smoothIntensity = DinoUtil.smoothMean(rawIntensity, 1)
		smoothIntensity = DinoUtil.smoothMean(
				DinoUtil.smoothMedian(rawIntensity, params.adv.hillSmoothMedianWindow),
				params.adv.hillSmoothMeanWindow
			)
		this
	}
	
	
	
	def smoothFullProfile:Seq[Double] = {
		val profile = new Array[Double](scanIndex.last - scanIndex.head + 1)
		for (i <- 0 until length)
			profile(scanIndex(i) - scanIndex.head) = smoothIntensity(i)
		for (i <- 0 until profile.length)
			if (profile(i) == 0) // should not need index guards since missing value is always surrounded by two values
				profile(i) = math.min(profile(i-1), profile(i+1)) - 1
			
		profile
	}
	
	
	
	def decompose(implicit params:DinosaurParams):List[Hill] = {
		if (params.adv.noHillSplit) return List(this)
		
		if (params.closeToDebug(501.1))
			{ val k = 1}
		
		splitIntoPair match {
			case Some((h1, h2)) =>
				h1.decompose ::: h2.decompose
			case None => 
				List(this)
		}
	}
	
	
	
	def splitIntoPair(implicit params:DinosaurParams):Option[(Hill, Hill)] = {
		if (params.adv.hillMaxMissing.value == 1) {
			val ints = smoothFullProfile
			val minPos = DinoUtil.localMinPositions(ints).sortBy(ints)
			for (pos <- minPos) {
				val leftMax = ints.take(pos).max
				val rightMax = ints.drop(pos+1).max
				val smallMax = math.min(leftMax, rightMax)
				if (smallMax / ints(pos) > params.adv.hillValleyFactor)
					return Some(splitAtScanIndex(scanIndex.head + pos))
			}
		} else {
			val ints = smoothIntensity
			val minPos = DinoUtil.localMinPositions(ints).sortBy(ints)
			for (pos <- minPos) {
				val leftMax = ints.take(pos).max
				val rightMax = ints.drop(pos+1).max
				val smallMax = math.min(leftMax, rightMax)
				if (smallMax / ints(pos) > params.adv.hillValleyFactor)
					return Some(splitAt(pos))
			}
		}
		
		return None
	}
	
	
	def splitAtScanIndex(ind:Int):(Hill, Hill) = 
		splitAt(scanIndex.indexWhere(_ >= ind))
	
	
	def splitAt(pos:Int):(Hill, Hill) = {
		val h0 = subHill(0, pos+1)
		h0.rightSplit = true
		val h1 = subHill(pos+1, length)
		h1.leftSplit = true
		(h0, h1)
	}
		
		
		
		
	def subHill(start:Int, end:Int):Hill = {
		val h = new Hill(
			centerMz.slice(start, end), 
			//minMz.slice(start, end),
			//maxMz.slice(start, end),
			rawIntensity.slice(start, end),
			scanIndex.slice(start, end))
		h.smoothIntensity = smoothIntensity.slice(start, end)
		h.maxMzWidth = maxMzWidth
		h
	}
}