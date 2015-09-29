package se.lth.immun

import se.jt.Params

import se.lth.immun.chem.IsotopeDistribution

class DinosaurAdvParams extends Params {

	import Params._
	
	// INTERNAL VALS
	// centroiding
	val subtractBackground 	= false 	## "(NOT IMPL) Remove spectral background"
	val backgroundQuantile 	= 0.0		## "(NOT IMPL) low abund quantile to remove"
	val intensityThreshold 	= 0.0		## "(NOT IMPL) hard intensity threshold below which peaks are removed"
	val maxIntensity		= false		## "Use max intensity for centroided peaks instead of the sum"
	val massEstPoints		= 3			## "Number of most abundant peaks to use for centroid determination"
	val massEstType:Estimation.EstType = Estimation.GaussianEst()
	
	
	// hill building
	val hillPPM = 8.0					## "Max ppm diff to allow in consecutive centroids for joining as a hill"
	val hillMaxMissing = 1				## "Max number of consecutive centroids that can be missing in a hill"
	val hillMinLength = 3				## "Min number of centroid in hill"
	//val hillMaxIntensity = false		## "Use max intensity for centroided peaks instead of the sum"
	val hillNBoots = 150				## "n bootstrap calculations for mz determination"
	val hillPeakFactor = 2				## "require hill smoothed endpoints to be < PeakFactor * apexIntensity"
	val hillPeakFactorMinLength = 40	## "only use hillPeakFactor on hills at least this long"
	val hillMzGuessLength = 3			## "number of recent mz values to use for live-guessing hill mz"
	val hillSmoothMedianWindow = 1		## "number of +-data points for the hill smooth median filter (applied before mean)"
	val hillSmoothMeanWindow = 1		## "number of +-data points for the hill smooth mean filter (applied after median)"
	val maxBootSize = 300				## "max number of samples to use in each bootstrap iteration during mz determination"
	val noHillSplit = false			## "Don't split hills based on intensity profile"
	val hillValleyFactor = 1.3 			## "Split hill if smallest surrounding local max is larger than valleyFactor * local min"
	val hillBatchSize = 100				## "Number of consecutive spectra to batch into non-parallel computation unit"
	
	// deisotoping
	val deisoCorr = 0.6					## "required cosine correlation between intensity profiles of hills in feature"
	val deisoSigmas = 5.0				## "number of m/z std-devs from expected to tolerate in feature isotopes"
	val deisoValleyFactor = 1.3			## "split feature if smallest surrounding local max is larger than valleyFactor * local min"
	val deisoOverlap = 3				## "minimum overlap required when using the averagine-overlap model"
	//val deisoCorrCalc:(Hill, Hill) => Double = Deisotoper.correlation(_,_)(this)
	val deisoCorrCalc:(Hill, Hill) => Double = _.corr(DinoUtil.cosineCorrelationFull)(_)
	val deisoDecompMaxSeeds = 100		## "max seeds checked for feature decomposition per round"
	val deisoAveragineStrategy:
		(Seq[Double], IsotopeDistribution, DinosaurParams) => Option[DinoUtil.Alignment] = 
			Cluster.averagineOverlap(_,_,_)
	val averagineCorr = 0.6				## "Required correlation with averagine"
	val averagineExplained = 0.5		## "Required probability mass of averagine explained"
	//val deisoAveragineStrategy:(Seq[Int], Seq[Double], se.lth.immun.chem.IsotopeDistribution, Int) => Deisotoper.IsotopePatternInds = Deisotoper.averagineFullCorr(_,_,_,_)
	
	val deisoBatchSize = 500			## "Number of hill clusters to batch into non-parallel computation unit"
	
	// chargePair
	val chargePairCorr = 0.6			## "required cosine correlation between intensity profiles of charge pairs"
	val massCalcNBoots = 150			## "n bootstrap calculations for mass calcuation"
	val chargePairPPM = 7.0				## "Max ppm diff to allow when assigning charge pairs"
	
	// report making
	val reportTargetIntMax = 2e3		## "Max value for spectral heatmaps unless bigger real max"
}