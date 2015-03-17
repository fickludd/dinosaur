package se.lth.immun

import akka.actor._
import java.io.File

import scala.collection.mutable.HashMap

object TopActor {
	case class AnalyzeFile(mzmlPath:File)
}

class TopActor(val params:DinosaurParams) extends Actor {

	import TopActor._
	import Centroider._
	import CentroidWorker._
	
	var ms1Specs:Option[Set[Int]] = None
	var centroidedSpecs = new HashMap[Int, Seq[CentroidPeak]]
	
	val centroider = context.actorOf(Props(new Centroider(params)), "centroider")
	
	def receive = {
		case AnalyzeFile(f) =>
			centroider ! Centroider.ImportAndCentroid(f)
		
		case FileParsed(specs, f) =>
			println("fininshed parsing file of "+specs.size+" spectra")
			ms1Specs = Some(specs)
			checkExit
			
		case CentroidingComplete(ind, res, t) =>
			centroidedSpecs += ind -> res
			params.centroidTime += t	
			checkExit
	}
	
	def checkExit = {
		if (centroidingFinished) {
			println("DONE!")
			if (params.profiling) {
				println("mzml parse time: "+Dinosaur.niceTiming(params.mzMLParseTime))
				println("  centroid time: "+Dinosaur.niceTiming(params.centroidTime))
				println(" avg centroid time: "+(params.centroidTime.toDouble / centroidedSpecs.size))
			}
			context.system.shutdown
		}
	}
	
	def centroidingFinished = 
		ms1Specs.nonEmpty && ms1Specs.get.size == centroidedSpecs.size
	
}