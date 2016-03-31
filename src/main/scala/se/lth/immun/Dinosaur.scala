package se.lth.immun

import se.jt.CLIApp
import java.util.Properties
import java.io.File
import java.io.Writer
import java.io.FileWriter
import java.io.BufferedWriter

import scala.util.Random

import se.lth.immun.xml.XmlReader
import se.lth.immun.mzml.MzML
import se.lth.immun.mzml.MzMLDataHandlers
import se.lth.immun.mzml.Spectrum

import se.lth.immun.mzml.ghost.GhostSpectrum


object Dinosaur extends CLIApp {

	def main(args:Array[String]):Unit = {
		
		var properties = new Properties
    	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
    	val name 		= properties.getProperty("pom.artifactId")
    	val version 	= properties.getProperty("pom.version")
    	val buildTime	= properties.getProperty("build.time")
    	
    	val params = new DinosaurParams(name, version)
    	
		params.startTime 	= System.currentTimeMillis
		
		failOnError(parseArgs(name, version, args, params, List("mzML"), None))
    	failOnError(parseParams(params.adv, params.advParams.value))
    	
		if (params.advHelp) {
			println(usage(name, version, args, params.adv, List("mzML"), None))
			System.exit(0)
		} 
		
		
		val (outDir, outName) = params.outBase
		
		println(name + " "+version + "    built:"+buildTime)
    	println("  mzML file: " + params.mzML.value)
    	println("    out dir: " + outDir)
    	println("   out name: " + outName)
		println()
		
		params.reportRandom = 
			if (params.reportSeed >= 0) 	new Random(params.reportSeed)
			else 							new Random
		
		val targets 	= (new TargetFile).read(params.targets, params.verbose)
		failOnError(params.setup(targets))
		
		val streamer 	= ReportStreamer(params, targets, outQcZipFile(params.outBase))
		
		val ff = new FeatureFinder(targets, streamer, params)
		val dinoResult = ff.analyzeMzML(new File(params.mzML))
		streamer.close
		
		println
		println("=== OUTPUT ===")
		var writeTime = 0L
		if (params.globalMode) {
			val w = initWriting(outFeatureFile(params.outBase), "global isotope pattern", params)
			writeTime += w.writeFeatureCsv(dinoResult, ff.reader.specTime)
		}
		
		if (params.targetMode) {
			val w = initWriting(outTargetFile(params.outBase), "targeted isotope pattern", params)
			writeTime += w.writeTargetCsv(dinoResult, ff.reader.specTime)
		}
		
		if (params.writeHills) {
			val w = initWriting(outHillFile(params.outBase), "hills", params)
			writeTime += w.writeHillCsv(dinoResult, ff.reader.specTime)
		}
		
		if (params.writeMsInspect) {
			val w = initWriting(outMsInspectFile(params.outBase), "msInspect", params)
			writeTime += w.writeMsInspectCsv(dinoResult, ff.reader.specTime)
		}
		
		if (params.writeBinary) {
			val binaryFile = outBinaryFile(params.outBase)
			println("  writing %s file %s".format("binary", binaryFile.toString))
			writeTime += OutWriter.writeBinary(binaryFile, dinoResult.patterns, ff.reader.specTime, params)
		}
		
		if (params.writeQuantML) {
			val mzqFile = outQuantMLFile(params.outBase)
			println("  writing %s file %s".format("mzQuantML", mzqFile.toString))
			writeTime += MzQuantML.write(mzqFile, dinoResult.patterns, ff.reader.specTime, params)
		}
		
		
		
		if (params.profiling) {
			println
			println("=== PROFILING ===")
			println("         mzml parse time: "+niceTiming(params.mzMLParseTime))
			println("           centroid time: "+niceTiming(params.centroidTime))
			println("        hill report time: "+niceTiming(params.hillReportTime))
			println("          deisotope time: "+niceTiming(params.deisotopeTime))
			println(" iso pattern report time: "+niceTiming(params.deisoReportTime))
			println("         mass calib time: "+niceTiming(params.massCalibTime))
			println("              write time: "+niceTiming(writeTime))
			println
			println("         deiso edge time: "+niceTiming(params.deisoEdgeTime))
			println("          deiso rip time: "+niceTiming(params.deisoRipTime))
			println("       deiso deconv time: "+niceTiming(params.deisoDeconvolveTime))
		}
		println
		println("total time: "+niceTiming(System.currentTimeMillis - params.startTime))
	}
	
	
	
	
	
	
	def toOutFile(ext:String)(base:(String, String)):File = 
		new File(base._1, base._2+"."+ext)
	
	def outFeatureFile = toOutFile("features.tsv") _
	def outTargetFile = toOutFile("targets.csv") _
	def outHillFile = toOutFile("hills.csv") _
	def outBinaryFile = toOutFile("features.bin") _
	def outQuantMLFile = toOutFile("mzq") _
	def outMsInspectFile = toOutFile("msInspect.tsv") _
	def outQcZipFile = toOutFile("qc.zip") _
	
	def initWriting(f:File, fileType:String, params:DinosaurParams) = {
		println("  writing %s file %s".format(fileType, f.toString))
		new OutWriter(params, new BufferedWriter(new FileWriter(f)))
	}
}