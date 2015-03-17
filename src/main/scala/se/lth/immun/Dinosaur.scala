package se.lth.immun

import se.jt.CLIApp
import java.util.Properties
import java.io.File

import akka.actor._

object Dinosaur extends CLIApp {

	
	def main(args:Array[String]):Unit = {
		
		var properties = new Properties
    	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
    	val name 		= properties.getProperty("pom.artifactId")
    	val version 	= properties.getProperty("pom.version")
    	
    	val params = new DinosaurParams(name, version)
    	
		params.startTime 	= System.currentTimeMillis
		
		failOnError(parseArgs(name, version, args, params, List("mzML"), None))
    	
		println(name + " "+version)
    	println("  mzML file: " + params.mzML.value)
		println()
		
		if (params.reportSeed >= 0)
			scala.util.Random.setSeed(params.reportSeed)
		
		if (params.nReport > 0) {
			val qcDir = new File("qc")
			if (!qcDir.exists)
				qcDir.mkdir
		}
		
		val system = ActorSystem("diana-system")
		val top = system.actorOf(Props(new TopActor(params)), "top-actor")
		top ! TopActor.AnalyzeFile(new File(params.mzML.value))
		
		system.awaitTermination
	}
}