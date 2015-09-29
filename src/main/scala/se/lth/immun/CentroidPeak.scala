package se.lth.immun

case class CentroidPeak(mz:Double, int:Double, minMz:Double, maxMz:Double) {
	def mzw = maxMz - minMz
}