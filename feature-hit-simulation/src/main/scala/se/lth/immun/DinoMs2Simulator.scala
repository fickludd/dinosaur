package se.lth.immun

import se.jt.CLIApp
import java.util.Properties

import java.io.File
import java.io.FileWriter
import java.io.BufferedWriter

import scala.io.Source
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.HashMap

import se.lth.immun.unimod.UniMod
import se.lth.immun.chem.Ion

object DinoMs2Simulator extends CLIApp {

	var properties = new Properties
	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
	val name 		= properties.getProperty("pom.artifactId")
	val version 	= properties.getProperty("pom.version")
	val buildTime	= properties.getProperty("build.time")
	
	val params = new DinoMs2Params
    
	case class MzZ(mz:Double, z:Int)
	trait Event{
		def mz:Double 
		def z:Int
		def rt:Double
	}
	case class ScrambleEvent(mz:Double, z:Int, rt:Double) extends Event
	
	case class Hit(
			file:String,
			spectrumID:Int,
			qValue:Double,
			peptide:String,
			intensity:Double,
			z:Int,
			rt:Double,
			precQuant:String,
			mz:Double
		) extends Event
		
	case class Feature(
			z:Int,
			mz:Double,
			file:String,
			rtStart:Double,
			rtApex:Double,
			rtEnd:Double,
			intensity:Double,
			peptide:String
		)
	
	def equals(e:Event)(f:Feature) = 
		e.z == f.z && 
		e.rt >= f.rtStart - params.rtTolerance && 
		e.rt < f.rtEnd + params.rtTolerance && 
		math.abs(e.mz - f.mz) < params.mzTol.absTolAt(e.mz)
		
	class FeaturesAndHits {
		var features = new ArrayBuffer[Feature]
		var hits = new ArrayBuffer[Hit]
		
		def sort = {
			features = features.sortBy(_.mz)
			hits = hits.sortBy(_.mz)
		}
	}
		
	case class Result(file:String, nHits:Int, nFeatures:Int, targetMatched:Double, scrambleMatched:Double) {
		def estimateFalseMatches:Double = {
			val k = (1-targetMatched) / (1-scrambleMatched)
			k * scrambleMatched
		}
		def estimateTrueMatches = targetMatched - estimateFalseMatches
	}
		
		
	
	def main(args:Array[String]):Unit = {
		
		failOnError(parseArgs(name, version, args, params, List("hits", "features"), None))
		
		print("reading hits...")
    	val hits = parseHits(params.hits)
    	println(" read "+hits.length)
    	print("reading features...")
    	val features = parseFeatures(params.features)
    	println(" read "+features.length)
    	
    	val grouped = groupByFile(features, hits)
    	println("grouped features and hits by file (%d groups)".format(grouped.size))
    	
    	val hitStatistics = 
    		for ((file, featAndHits) <- grouped) yield 
    			Result(
    				file, 
    				featAndHits.hits.length,
    				featAndHits.features.length,
    				matchFeaturesToEvent(featAndHits.features, featAndHits.hits),
					matchFeaturesToEvent(featAndHits.features, scramble(featAndHits.hits))
				)
    	
		println("avg. est true matches: %.5f".format(hitStatistics.map(_.estimateTrueMatches).sum / hitStatistics.size))
		
    	printResults(hitStatistics)
	}
	
	
	
	
	
	def printResults(res:Iterable[Result]) = {
		val w = new BufferedWriter(new FileWriter(params.outFile))
		
		def row(xs:Any*) =
			w.write(xs.mkString("\t") + "\n")
		
		row("hitsFile", "featureFile", "file", "nHits", "nFeatures", "mzTolerance", "rtTolerance", "targetMatches", "scrambledMatches", "estFalseMatches", "estTrueMatches")
		for (r <- res) 
			row(params.hits.value, params.features.value, r.file, r.nHits, r.nFeatures, params.mzTol, params.rtTolerance.value, r.targetMatched, r.scrambleMatched, r.estimateFalseMatches, r.estimateTrueMatches)
			
		w.close
	}
	
	
	def scramble(hits:Seq[Hit]):Seq[Event] = {
		val mzzs 	= params.random.shuffle(hits.map(h => MzZ(h.mz, h.z)))
		val rts 	= params.random.shuffle(hits.map(_.rt))
		
		val scrambled = 
			for (i <- 0 until rts.length) yield ScrambleEvent(mzzs(i).mz, mzzs(i).z, rts(i))
		
		scrambled.sortBy(_.mz)
	}
	
	
	def matchFeaturesToEvent(features:Seq[Feature], events:Seq[Event]):Double = {
		var count = 0
		val t0 = System.currentTimeMillis
		
		if (params.matchProcedure.value == "slow")
			count = events.count(e => features.exists(equals(e)))
		else { 
			var fiLow = 0
			var fiHigh = 0
			
			for (e <- events) {
				val mzTol = params.mzTol.absTolAt(e.mz)
				while (fiLow < features.length && features(fiLow).mz < e.mz - mzTol) fiLow += 1
				fiHigh = math.max(fiHigh, fiLow)
				while (fiHigh < features.length && features(fiHigh).mz < e.mz + mzTol) fiHigh += 1
				
				var fi = fiLow
				var matched = false
				while (fi < fiHigh) {
					if (equals(e)(features(fi))) {
						matched = true
						fi = fiHigh
					} else
						fi += 1
				}
				
				if (matched) count += 1
				if (params.debugMode)
					e match {
						case s:ScrambleEvent => {}
						case h:Hit =>
							if (h.precQuant == "null" && matched)
								println("Matched %s even though Proteios didn't".format(h))
							else if (h.precQuant != "null" && !matched)
								println("Proteios matched %s by I didn't".format(h))
					}
			}
		}
		
		val t1 = System.currentTimeMillis
		if (params.debugMode)
			println("matching complete in %10d ms".format(t1 - t0))
		 
		count.toDouble / events.length
	}
	
	
	def groupByFile(features:Seq[Feature], hits:Seq[Hit]) = {
		val map = new HashMap[String, FeaturesAndHits]
		for (f <- features) {
			if (!map.contains(f.file))
				map += f.file -> new FeaturesAndHits
			map(f.file).features += f
		}
		for (h <- hits) {
			if (!map.contains(h.file))
				map += h.file -> new FeaturesAndHits
			map(h.file).hits += h
		}
		for (featAndHits <- map.values) featAndHits.sort 
		map
	}
	
	
	def parseHits(path:String) = {
		(for (line <- Source.fromFile(new File(path)).getLines.drop(1)) yield {
			val p = line.split("\t").map(_.trim)
			
			val file = p(1).split("=").last.dropRight(6)
			val spectrumID = p(2).toInt
			val qValue = p(3).toDouble
			val peptide = p(6)
			val intensity = p(11).toDouble
			val z = p(12).toInt
			val rt = p(13).toDouble
			val precQuant = p(14)
			
			val molecule = UniMod.parseUniModSequence(fixPeptide(peptide))
			
			Hit(file, spectrumID, qValue, peptide, intensity, z, rt, precQuant, Ion.mz(molecule, z))
		}).toSeq
	}
	
	def fixPeptide(x:String):String = {
		val rawPep = x.takeWhile(_ != ' ')
		val ptmPep = 
			if (x.contains("Oxidation (M) @")) {
				val mods = x.dropWhile(_ != '@').split("delta:").head.zipWithIndex.filter(_._1 == '1').map(t => ("(UniMod:35)", t._2 - 2))
				mergeRawPepWithMods(rawPep.zipWithIndex.toList, mods.toList)
			} else if (x.contains("delta:")) {
				val modStr = x.dropWhile(_ != ' ').split("delta:").head.trim
				val mods = modStr.split(" ").map(parseMod).sortBy(_._2)
				mergeRawPepWithMods(rawPep.zipWithIndex.toList, mods.toList)
			} else
				rawPep
		ptmPep.replaceAll("C", "C(UniMod:4)")
	}
	
	def mergeRawPepWithMods(rawPep:List[(Char, Int)], mods:List[(String, Int)]):String = {
		rawPep match {
			case Nil =>
				mods.map(_._1).mkString
			case (aa, aaPos)::aas =>
				mods match {
					case Nil =>
						rawPep.map(_._1).mkString
					case (mod, mPos)::more =>
						if (mPos <= aaPos)
							mod + mergeRawPepWithMods(rawPep, more)
						else if (mPos == aaPos) 
							aa + mod + mergeRawPepWithMods(aas, more)
						else
							aa + mergeRawPepWithMods(aas, mods)
				}
		}
	}
	
	def parseMod(modStr:String):(String, Int) = {
		val atIndex = modStr.indexOf('@')
		val pos = modStr.substring(atIndex+2).toInt
		val mod = modStr.take(atIndex+2)
		try {
			if (mod.startsWith("42.01057@") && pos == 1)
				("(UniMod:1)", 0)
			else
				mod match {
					case "-18.01056@E" if pos == 1 => ("(UniMod:27)", 0)
					case "-17.02655@Q" if pos == 1 => ("(UniMod:28)", 0)
					case "15.99492@M" => ("(UniMod:35)", pos)
					case "-17.02655@C" => ("(UniMod:28)", 0)
				}
		} catch {
			case e:Exception =>
				println("error parsing '%s'".format(modStr))
				throw e
		}
	}
	
	def parseFeatures(path:String) = {
		(for (line <- Source.fromFile(new File(path)).getLines.drop(1)) yield {
			val p = line.split("\t").map(_.trim)
			
			val z = p(0).toInt
			val mz = p(1).toDouble
			val file = p(2).split("=").last.dropRight(6)
			val rtStart = p(5).toDouble
			val rtApex = p(4).toDouble
			val rtEnd = p(6).toDouble
			val intensity = p(7).toDouble
			val peptide = p(8)
			
			Feature(z, mz, file, rtStart, rtApex, rtEnd, intensity, peptide)
		}).toSeq
	}
}