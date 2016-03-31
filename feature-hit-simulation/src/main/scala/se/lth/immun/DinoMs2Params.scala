package se.lth.immun

import se.jt.Params
import scala.util.Random
import java.io.File

class DinoMs2Params extends Params {

	import Params._
	
	trait Tolerance { def absTolAt(mz:Double):Double }
	case class AbsTol(dx:Double) extends Tolerance { 
		def absTolAt(mz:Double) = dx
		override def toString = "%.4f Da".format(dx)
	}
	case class PPMTol(ppm:Double) extends Tolerance {
		def absTolAt(mz:Double) = mz / 1e6 * ppm
		override def toString = "%.1f ppm".format(ppm)
	}
	
	val hits = ReqString("Proteios hits export file describing all MS/MS IDs")
	val features = ReqString("Proteios features export file describing all used features IDs")
	
	val mzTolerance = "10 ppm"			## "Tolerance in m/z to accept feature-hit match (ppm or Da)"
	val rtTolerance = 0.0				## "Tolerance in minutes outside feature rt start and end"
	val seed = 1					## "Seed for hits scrambling"
	val matchProcedure = "slow"		## "Procedure for match calculations (slow | fast)"
	val epsilon = 0.0001
	val debugMode = false			## "Set to compare Proteios and current matching"
	
	lazy val random = new Random(seed)
	lazy val mzTol:Tolerance = {
		val l = mzTolerance.toLowerCase
		if (l.endsWith("ppm")) 
			PPMTol(l.dropRight(3).trim().toDouble)
		else if (l.endsWith("da")) 
			AbsTol(l.dropRight(2).trim().toDouble)
		else
			throw new Exception("Unable to parse mzTolerance '%s'".format(mzTolerance.value))
	}
	
	def outFile = 
		if (hits.value.toLowerCase.endsWith(".tsv"))
			new File(hits.value.dropRight(4) + ".matched.tsv")
		else
			new File(hits.value + ".matched.tsv")
}