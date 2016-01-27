package se.lth.immun

import java.io.File
import akka.actor._
import se.jt.CLIBar
import scala.util.Random
import scala.util.Try
import scala.concurrent.Await
import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global
import scala.concurrent.duration.Duration
import scala.concurrent.duration.TimeUnit
import akka.pattern.ask
import akka.util.Timeout
import scala.collection.mutable.HashSet
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.Queue
import java.util.concurrent.TimeUnit

case class ImTheCustomer()
case class DinosaurResult(
		patterns:Seq[IsotopePattern], 
		targetMatches:Seq[TargetMatcher.Match], 
		massCalib:MzIntensity => Double
	)

class FeatureFinder( 
		targets:Seq[Target],
		streamer:ReportStreamer,
		implicit val params:DinosaurParams
) extends Timeable {

	import Deisotoper._
	
	val actorSystem 	= ActorSystem("actor-system")
	val actorInbox 		= Inbox.create(actorSystem)
	
	val reader				= new MzMLReader(params, streamer)
	val hillBuilderActor 	= actorSystem.actorOf(Props(new HillBuilderActorParallel(params)))
	val deisotoper 			= actorSystem.actorOf(Props(new Deisotoper(params)))
	val targetDeisotoper	= new TargetDeisotoper(params)
	val chargePairer 		= new ChargePairer(params)
	val nonLinMassCalibration = new NonLinMassCalibration(params)
	
	
	
	def await[T](pf:PartialFunction[AnyRef, T]):T = {
		var res:Option[T] = None
		while (res.isEmpty)
			try {
				actorInbox.receive(Duration.create(1, TimeUnit.SECONDS )) match {
					case str:String => if (params.verbose) print(str)
					case a:AnyRef => 
						if (pf.isDefinedAt(a)) 
							res = Some(pf(a))
				}
			} catch {
				case e:java.util.concurrent.TimeoutException => {}
			}
		res.get
	}
	
	def awaitHillBuilding = await({ case BuiltHills(hills) => hills })
	def awaitDeisotoping = await({ case DeisotopingComplete(clusters, patterns) => (clusters, patterns) })
	
	
	
	
	def analyzeMzML(f:File):DinosaurResult = {
		timer.start
		
		// READ SPECTRA, CENTROID THEM AND BUILD HILLS
		reader.read(f, (a:AnyRef) => actorInbox.send(hillBuilderActor, a))
		//actorInbox.send(hillBuilderActor, Ms1Count(reader.ms1Index))
		val completedHills = awaitHillBuilding 
		
		params.mzMLParseTime = timer.click
		
		val allHills = completedHills.sortBy(_.scanIndex.last).toArray
		println("all hills, n="+allHills.length)
		println("hill checkSum = "+allHills.zipWithIndex.map(t => t._1.total.centerMz.toLong * t._2).sum)
		val hills = allHills.filter(h => 
						h.length < params.adv.hillPeakFactorMinLength || 
						(h.smoothIntensity.head * params.adv.hillPeakFactor < h.apex.intensity && 
						 h.smoothIntensity.last * params.adv.hillPeakFactor < h.apex.intensity)
					)
		
		println("peaky hills, n="+hills.length)
		println("peaky hill checkSum = "+hills.zipWithIndex.map(t => t._1.total.centerMz.toLong * t._2).sum)
		printHillHistogram(hills)
		
		println("writing hill reports...")
		HillReport.createReports(hills, reader.ghostSpecMap, streamer, params)
		println("hill reports written")
		
		params.hillReportTime = timer.click
		
		val (clusters, isotopes) = 
			if (params.globalMode) {
				actorInbox.send(deisotoper, Deisotope(hills, reader.specTime))
				awaitDeisotoping
			} else 
				targetDeisotoper.deisotope(hills, targets, reader.specTime)
		println("deisotoping complete")
		println("isotopes, n="+isotopes.length)
		
		params.deisotopeTime = timer.click
		
		if (params.globalMode) {
			println("writing isotope pattern reports...")
			DeisotopeReport.createReports(isotopes, reader, streamer)
			IsotopeDistributionReport.createReports(isotopes, streamer)
			println("isotope pattern reports written")
		}
		
		params.deisoReportTime = timer.click
		
		val pairs = chargePairer.pairPatterns(isotopes)
		val massCalib = 
			Try(nonLinMassCalibration.findRecalibrationFunction(pairs))
		//NonLinMassCalibReport.createReports(isotopes, massCalib, streamer)
		
		params.massCalibTime = timer.click
		
		actorSystem.shutdown
		
		val targetMatches =
			if (targets.nonEmpty) TargetMatcher.matchTargets(targets, isotopes, reader.specTime)
			else Nil
		
		if (params.targetsToReport) 
			TargetReport.createReports(targetMatches, isotopes, hills, reader, streamer)
		
		DinosaurResult(
				isotopes, 
				targetMatches, 
				massCalib.getOrElse((mzInt:MzIntensity) => mzInt.mz)
			)
	}
	
	
	def printHillHistogram(hills:Seq[Hill]) = {
		
		def hist(str:String, f:Int => Boolean) = 
			println("%10s %8d".format(str, hills.count(h => f(h.length))))
		
		def countVal(n:Int) =
			hist(""+n, _ == n)
			
		def countRange(a:Int, b:Int) =
			hist(a+"-"+b, l => l >= a && l < b)
			
		println("  nScans    nHills")
		countVal(2)
		countVal(3)
		countVal(4)
		countRange(5, 10)
		countRange(10, 20)
		countRange(20, 50)
		countRange(50, 100)
		countRange(100, 200)
		countRange(200, 500)
		countRange(500, 1000)
		countRange(1000, 2000)
		countRange(2000, 5000)
		countRange(5000, 10000)
		hist(">10000", l => l > 10000)
	}
}