package se.lth.immun

import scala.collection.mutable.ArrayBuffer
import java.io.File
import scala.io.Source


case class Target(
		val mz:Double, 
		val z:Int, 
		val mzDiff:Double, 
		val rtStart:Double, 
		val rtEnd:Double, 
		val minApexInt:Double, 
		val id:String,
		val group:Option[String]
	)

class TargetFile {
	
	val errs = new ArrayBuffer[String]
	val columns = Array("mz", "charge", "mzDiff", "rtStart", "rtEnd", "minApexInt", "id", "group")
	
	def indexOfCol(header:Seq[String], col:String) = {
		val i = header.indexOf(col)
		if (i < 0)
			errs += "Could not find column '%s' in feature file".format(col)
		i
	}
	
	
	
	def read(path:String, verbose:Boolean):Seq[Target] = {
		val targets = read(path)	
		if (verbose) {
			println("Target Reports (n=%d):".format(targets.length))
			println("   mz     rtStart  rtEnd   minApexInt   id")
			for (tr <- targets) 
				println("%8.4f %6.2f %6.2f %6.1e %s %s".format(tr.mz, tr.rtStart, tr.rtEnd, tr.minApexInt, tr.group.getOrElse(""), tr.id))
		}
		targets
	}
	
	
	
	def read(path:String):Seq[Target] = {
		if (path != "") {
			val f = new File(path)
			if (f.exists)
				parse(f)
			else 
				throw new Exception("Could not read report targets from file '%s'. Needed columns are: %s".format(path, columns.mkString(",")))
		} else Nil
	}
	
	
	def parse(f:File) = {
		var headerParsed = false
		
		var iMZ = -1
		var iZ = -1
		var iMZ_DIFF = -1
		var iRT_START = -1
		var iRT_END = -1
		var iMIN_INT = -1
		var iID = -1
		var iGROUP = -1
		
		var res = new ArrayBuffer[Target]
		try {
			for (line <- Source.fromFile(f).getLines) {
				val cols = line.split("\t")map(_.trim)
				if (!headerParsed) {
					val header = cols
					iMZ = indexOfCol(header, columns(0))
					iZ = indexOfCol(header, columns(1))
					iMZ_DIFF = indexOfCol(header, columns(2))
					iRT_START = indexOfCol(header, columns(3))
					iRT_END = indexOfCol(header, columns(4))
					iMIN_INT = indexOfCol(header, columns(5))
					iID = indexOfCol(header, columns(6))
					iGROUP = indexOfCol(header, columns(7))
					headerParsed = true
				} else {
					res += Target(
						cols(iMZ).toDouble,
						cols(iZ).toInt,
						cols(iMZ_DIFF).toDouble,
						cols(iRT_START).toDouble,
						cols(iRT_END).toDouble,
						cols(iMIN_INT).toDouble,
						cols(iID),
						if (iGROUP >= 0) Some(cols(iGROUP)) else None
					)
				}
			}
		} catch {
			case e:Exception => {e.printStackTrace()}
		}
		res
	}
}