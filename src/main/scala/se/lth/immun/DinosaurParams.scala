package se.lth.immun

import se.jt.Params

class DinosaurParams(val name:String, val version:String) extends Params {

	import Params._
	
	// USER EXPOSED PARAMS
	val verbose 		= false 	## "increase details in output"
	val concurrency 	= 1 		## "the number of assays to analyze in parallel"
	val profiling 		= false 	## "set to enable CPU profiling"
	val force			= false		## "ignore missing mzML params"
	val nReport			= 10		## "number of random assay to export control figure for"
	val reportSeed		= -1L		## "seed to use for report assay selection (<0 means random)"
	
	val mzML = ReqString("The shotgun MzML file to analyze")
	
	
	// INTERNAL VALS
	// centroiding
	val subtractBackground = false
	val backgroundQuantile = 0
	val intensityThreshold = 0
	val maxIntensity		= false
	val massEstPoints		= 3
	val massEstType:CentroidWorker.MassEstType = CentroidWorker.GaussianMassEst()
	
	// INTERNAL VARS
	var startTime		= 0L
	var mzMLParseTime	= 0L
	var centroidTime 	= 0L
	
}