package se.lth.immun

import akka.actor._
import java.io.File
import java.io.FileReader
import java.io.FileInputStream
import java.io.BufferedReader
import java.io.InputStreamReader
import java.util.zip.GZIPInputStream

import se.lth.immun.xml.XmlReader
import se.lth.immun.mzml.MzML
import se.lth.immun.mzml.MzMLDataHandlers
import se.lth.immun.mzml.Spectrum

import se.lth.immun.mzml.ghost.GhostSpectrum

import scala.util.Random
import scala.collection.mutable.HashSet
import scala.collection.mutable.HashMap

object Centroider {
	case class ImportAndCentroid(f:File)
	case class FileParsed(ms1Specs:Set[Int], f:File)
}

class Centroider(val params:DinosaurParams) extends Actor {

	import Centroider._
	import CentroidWorker._
	
	val specs = new HashSet[Int]
	var reportSpectra:Seq[Int] = _
	var customer:ActorRef = _ 
	
	def receive = {
		case ImportAndCentroid(f) =>
			customer = sender
			
			val dh = new MzMLDataHandlers(
				setupDataStructures,
				handleSpectrum,
				nc => {},
				c => {})
			
			val xr = new XmlReader(
					if (f.getName.toLowerCase.endsWith(".gz"))
						new BufferedReader(new InputStreamReader(
							new GZIPInputStream(new FileInputStream(f))))
					else
						new BufferedReader(new FileReader(f))
				)
			
			xr.force = params.force
			val mzML = MzML.fromFile(xr, dh)
			sender ! FileParsed(specs.toSet, f)
			params.mzMLParseTime = System.currentTimeMillis - params.startTime
			
			
	}
	
	def setupDataStructures(numSpec:Int) = {
		val specs = (0 until numSpec).toTraversable
		reportSpectra = Random.shuffle(specs).take(params.nReport).toSeq.sorted
	}
	
	def handleSpectrum(s:Spectrum):Unit = {
		val gs = GhostSpectrum.fromSpectrum(s)
		
		if (gs.msLevel == 1) {
			specs += gs.spectrum.index
			
			var reportCurrent = false
			if (reportSpectra.nonEmpty && gs.spectrum.index > reportSpectra.head) {
				reportCurrent = true
				reportSpectra = reportSpectra.tail
			}
			
			val c = context.actorOf(Props(new CentroidWorker(params)))
			c ! CentroidWorker.CentroidSpectrum(gs, customer, reportCurrent)
		}
	}
}