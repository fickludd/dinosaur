package se.lth.immun

import se.jt.Params
import scala.util.Random
import java.io.File

import scala.collection.mutable.ArrayBuffer

class DinosaurParams(val name:String, val version:String) extends Params {

	import Params._
	
	// USER EXPOSED PARAMS
	val verbose 		= false 	## "increase details in output"
	val concurrency 	= 2 		## "the number of assays to analyze in parallel"
	val profiling 		= false 	## "set to enable CPU profiling"
	val force			= false		## "ignore missing mzML params"
	
	val writeHills		= false		## "set to output csv file with all hills assigned to isotope patterns"
	val writeMsInspect	= false		## "set to output MsInspect feature csv file"
	val writeBinary		= false		## "set to output binary MSFeatureProtocol file"
	val writeQuantML	= false		## "set to output mzQuantML file" 
	val outDir			= ""		## "output directory (by default same as input mzML)"
	val outName			= ""		## "basename for output files (by default same as input mzML)"
	
	val nReport			= 10		## "number of random assay to export control figure for"
	val reportSeed		= -1L		## "seed to use for report assay selection (<0 means random)"
	val reportDeisoMzHeight	= 15.0	## "mz range in deisotoper reports"
	val reportHighRes	= false		## "generate high-resolution plot trail when supported (for print)"
	val reportTargets	= false		## "set to create a special report figure for each target"
	val targets = ""				## "path to isotope patterns target file (not used by default)"
	val targetPreference = "rt"		## "if multiple isotope patterns fit target, take the closest rt apex (rt) or the most intense (intensity)"
	val mode = "global"				## "analysis mode: global or target. Global mode reports all isotope patterns, targeted only those matching targets."
	val zipQcFolder = false			## "set to zip the entire qc folder on algorithm completion"
	
	val minCharge = 1				## "min searched ion charge"
	val maxCharge = 6				## "max searched ion charge"
	
	val advParams = ""				## "path to adv param file"
	val advHelp = false				## "set to output adv param file help and quit"
	val adv = new DinosaurAdvParams
	
	val mzML = ReqString("The shotgun MzML file to analyze")
	
	def outBase = {
		val mzMLFile = new File(mzML)
		val dir = 
			if (outDir.value != "") outDir.value
			else mzMLFile.getParent
		val name =
			if (outName.value != "") outName.value
			else stripExt(mzMLFile.getName)
		(dir, name) 
	}
	
	val ALLOWED_MODES = Array("global", "target")
	def globalMode = mode.value == "global"
	def targetMode = mode.value == "target"
	
	val ALLOWED_TARGET_PREF = Array("rt", "intensity")
	
	def setup(targets:Seq[Target]) = {
		if (targets.nonEmpty && (globalMode || reportTargets)) { 
			spectrumBacklogSize = 1000000
			targetsToReport = true
		}
		
		val errs = new ArrayBuffer[String]
		if (ALLOWED_MODES.forall(_ != mode.value))
			errs += "Unsupported mode '%s'. Should be one of: %s".format(mode.value, ALLOWED_MODES.mkString(", "))
			
		if (ALLOWED_TARGET_PREF.forall(_ != targetPreference.value))
			errs += "Unsupported target preference '%s'. Should be one of: %s".format(mode.value, ALLOWED_TARGET_PREF.mkString(", "))
		
		//println("SPECTRUM BACKLOG: "+spectrumBacklogSize)
		errs
	}
		
	def stripExt(path:String) =
		if (path.toLowerCase.endsWith(".mzml"))
			path.dropRight(5)
		else if (path.toLowerCase.endsWith(".mzml.gz"))
			path.dropRight(8)
		else
			path
	
	
	// output write
	val outQuote = "\""
	val outSep = "\t"
		
	// verbose
	val freqHillRefine = 10000
	val freqFindEdge = 10000
	val freqClusterDeconvolve = 1000
	var spectrumBacklogSize = 150
	var targetsToReport = false
	
	// INTERNAL VARS
	var startTime		= 0L
	var mzMLParseTime	= 0L
	var centroidTime 	= 0L
	var hillReportTime = 0L
	var deisotopeTime 	= 0L
	var deisoReportTime = 0L
	var massCalibTime 	= 0L
	var deisoEdgeTime 	= 0L
	var deisoRipTime	= 0L
	var deisoDeconvolveTime = 0L
	
	var reportRandom:Random = _
	
	// DEBUG help
	var debugMz = -1.0
	val debugMzDiff = 0.01
	var debugRtMinStart = -1.0
	var debugRtMinEnd = 10000000.0
	
	var debugHill:Option[Hill] = None
	var debugHillInd:Option[Int] = None
	var debugCluster:Option[Int] = None
	
	//setDebug(506.78, 45.1, 45.7)
	def setDebug(mz:Double, rtStart:Double, rtEnd:Double) = {
		debugMz = mz
		debugRtMinStart = rtStart
		debugRtMinEnd = rtEnd
	}
	
	def findDebugHill(hills:Seq[Hill], specTime:Seq[Double]) = {
		debugHill = hills.find(h => 
			specTime(h.apex.scanIndex) > debugRtMinStart && specTime(h.apex.scanIndex) < debugRtMinEnd && 
			closeToDebug(h.total.centerMz))
		debugHillInd = debugHill.map(h => hills.indexOf(h))
	}
	
	def findDebugCluster(clusters:Seq[Seq[Cluster.Edge]]) =
		debugCluster = debugHillInd.map(i => clusters.indexWhere(_.exists(e => e.i == i || e.j == i)))
	
	def closeToDebug(mz:Double) =
		math.abs(mz - debugMz) < debugMzDiff
		
	def inDebugTimeWindow(rt:Double) =
		rt > debugRtMinStart - 1 && rt < debugRtMinEnd + 1
		
	def atDebugTime(rt:Double) =
		rt > debugRtMinStart && rt < debugRtMinEnd
}
