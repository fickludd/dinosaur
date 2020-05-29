package se.lth.immun

import se.lth.immun.chem.Constants

case class IsotopePattern(
		val hills:Seq[Hill], 
		val isotopeOffset:Int, 
		val mostAbundNbr:Int, 
		val z:Int, 
		val averagineCorr:Double,
		val bootstrap:SplittableBootstrap
)(implicit params:DinosaurParams) {
	
	def length 		= hills.length
	val startInd 	= hills.map(_.scanIndex.head).min
	val endInd 		= hills.map(_.scanIndex.last).max
	lazy val totalProfile 	= hillsProfile(hills)
	
	val (mass, massError) = bootstrapMass(params.adv.massCalcNBoots)
	
	lazy val apexHill	= hills.maxBy(_.apex.intensity)
	lazy val intensity = hills.map(_.total.intensity).sum
	def startRt(specTime:Seq[Double]) 	= hills.map(_.startRt(specTime)).min
	def endRt(specTime:Seq[Double]) 	= hills.map(_.endRt(specTime)).max
	
	
	
	def bootstrapMass(nBoot:Int) = {
		val rawIntensities = hills.flatMap(_.rawIntensity)
		val centerMz = hills.flatMap(_.centerMz)
		val isotopeNbrs = hills.zipWithIndex.flatMap(t => t._1.centerMz.map(_ => t._2+isotopeOffset))
		val n = hills.map(_.length).sum
		
		def isoMass(mz:Double, isotopeNbr:Int) =
			z*(mz - Constants.PROTON_WEIGHT) - DinoUtil.ISOTOPE_PATTERN_DIFF * isotopeNbr
	  
		val mPreCalc = (0 until rawIntensities.length).map(i => (rawIntensities(i) * isoMass(centerMz(i), isotopeNbrs(i)) ) )
		val ms = 
			for (i <- 0 until nBoot) 
			yield averageMass(
					rawIntensities,
					bootstrap.set(math.min(params.adv.maxBootSize, n), n, nBoot),
					mPreCalc
				)
		
		val mass = ms.sum / nBoot
		val massError = math.sqrt(ms.map(m => (mass - m)*(mass - m)).sum / nBoot)
		(mass, massError)
	}
	
	
	
	def averageMass(
			rawInts:Seq[Double], 
			indices:Seq[Int],
			mPreCalc:Seq[Double]
	):Double = {
	  var m = 0.0
	  var norm = 0.0
	  var i = 0
	  while (i < indices.length) {
		  m += mPreCalc(indices(i))
		  norm += rawInts(indices(i))
		  i += 1
		}
		m / norm
	}
	
	
	
	def hillsProfile(hills:Seq[Hill]):HillProfile = {
		val inds = hills.map(_.scanIndex.head).min to hills.map(_.scanIndex.last).max
		val ints = new Array[Double](inds.length)
		for {
			h <- hills
			(ind, int) <- h.scanIndex.zip(h.smoothIntensity)
		} (
			ints(inds.indexOf(ind)) += int
		)
		HillProfile(inds.head, ints)
	}
}
