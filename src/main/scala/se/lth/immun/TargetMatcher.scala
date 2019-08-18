package se.lth.immun

object TargetMatcher {
	case class Match(target:Target, pattern:Option[IsotopePattern])
	
	def matchTargets(
			targets:Seq[Target], 
			patterns:Seq[IsotopePattern],
			specTime:Seq[Double]
	)(implicit params:DinosaurParams) = {
	  val sortedPatterns = patterns.toArray.sortBy(_.hills.head.total.centerMz)
	  val patternMzs = sortedPatterns.map(_.hills.head.total.centerMz)
	  val accurateApexRts = sortedPatterns.map(_.apexHill.accurateApexRt(specTime))
	  var i = 0
		for (t <- targets) yield {  	  
  	  if (i % 10000 == 0 && params.verbose) {
		    println("Matching target " + i)
		  }
		  i += 1
		  
  	  val (patternsStartIndx, patternsEndIndx) = DinoUtil.getMinMaxIndx(patternMzs, t.mz, t.mzDiff)		  
			val candidates = for {
			  i <- patternsStartIndx until patternsEndIndx
			  val ip = sortedPatterns(i)
			  val rt = accurateApexRts(i)
			  if (rt > t.rtStart && rt < t.rtEnd && 
			        ip.apexHill.apex.intensity >= t.minApexInt && t.z == ip.z)
		  } yield (ip, rt)
		  
			val targetMidRt = (t.rtStart + t.rtEnd) / 2
			val best =
				if (candidates.nonEmpty) 
					Some(
						if (params.targetPreference.value == "close")
							candidates.minBy(p => math.abs(p._2 - targetMidRt))._1
						else 
							candidates.maxBy(_._1.apexHill.apex.intensity)._1
					)
				else None
			Match(t, best)
		}
  }
}
