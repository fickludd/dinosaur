package se.lth.immun

import java.io.File
import java.io.FileReader
import java.io.FileInputStream
import java.io.BufferedReader
import java.io.InputStreamReader
import java.util.zip.GZIPInputStream

import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer

import se.jt.CLIBar

import se.lth.immun.xml.XmlReader
import se.lth.immun.mzml.MzML
import se.lth.immun.mzml.MzMLDataHandlers
import se.lth.immun.mzml.Spectrum
import se.lth.immun.mzml.ghost.GhostSpectrum
import se.lth.immun.mzml.ghost.GhostException

class MzMLReader(params:DinosaurParams, streamer:ReportStreamer) {
	
	
	val centroider		= new Centroider(params, streamer)
	val readBar 		= new CLIBar
	val ghostSpecMap 	= new HashMap[Int, Array[GhostSpectrum]]
	val specTime 		= new ArrayBuffer[Double]
	
	var reportSpectra:Seq[Int] = _
	var send:AnyRef => Unit = a => {}
	
	var nSpec = 0
	var ms1Index = 0
	var specBackLog = List[GhostSpectrum]()
	
	
	val dh = new MzMLDataHandlers(
			setupDataStructures,
			handleSpectrum,
			nc => {},
			c => {})
	
	
	
	def ms1IndexAtRt(rt:Double) = 
		specTime.indexOf(specTime.minBy(ms1rt => math.abs(rt - ms1rt)))
	
		
	
	def indexOverlapRt(ms1IndexLow:Int, ms1IndexHigh:Int, rtStart:Double, rtEnd:Double) =
		specTime(ms1IndexHigh) >= rtStart && specTime(ms1IndexLow) < rtEnd
		
	
	
	def read(f:File, send:AnyRef => Unit) = {
		val xr = new XmlReader(
				if (f.getName.toLowerCase.endsWith(".gz"))
					new BufferedReader(new InputStreamReader(
						new GZIPInputStream(new FileInputStream(f))))
				else
					new BufferedReader(new FileReader(f))
			)
		
		xr.force = params.force
		this.send = send
		MzML.fromFile(xr, dh)
		
		if (!params.verbose) 
			print(readBar.update(nSpec, nSpec))
			
		send(Ms1Count(ms1Index))
	}
	
	
	
	def setupDataStructures(numSpec:Int) = {
		if (!params.verbose) {
			println(readBar.reference)
			print(readBar.update(0, numSpec))
		} 
		
		nSpec = numSpec
		val specs = (0 until numSpec).toTraversable
		reportSpectra = params.reportRandom.shuffle(specs).take(params.nReport).toSeq.sorted
		
		send(BuildFromSpectra(numSpec))
	}
	
	
	
	def handleSpectrum(s:Spectrum):Unit = {
		val gs = 
			try {
				GhostSpectrum.fromSpectrum(s)
			} catch {
				case ge:GhostException =>
					println("Error parsing spectrum '%d'".format(s.index))
					throw ge
			}
			
		if (gs.intensities.isEmpty)
			return
		
		if (!params.inDebugTimeWindow(gs.scanStartTime))
			return
			
		if (gs.msLevel == 1) {
			specTime += gs.scanStartTime
		
			specBackLog = gs :: specBackLog.take(params.spectrumBacklogSize - 1)
			
			var reportCurrent = false
			if (reportSpectra.nonEmpty && gs.spectrum.index > reportSpectra.head) {
				reportCurrent = true
				reportSpectra = reportSpectra.tail
				ghostSpecMap(ms1Index - specBackLog.length + 1) = specBackLog.reverse.toArray
			}
			
			val (cspec, t) = centroider.centroidSpectrum(gs, reportCurrent)
			send(AddSpectrum(gs.spectrum.index, ms1Index, cspec))
			ms1Index += 1
		}
			
		if (!params.verbose) 
			print(readBar.update(math.max(0, gs.spectrum.index - nSpec), nSpec))
	}
}