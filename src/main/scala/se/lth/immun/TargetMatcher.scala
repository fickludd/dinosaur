package se.lth.immun

object TargetMatcher {
	case class Match(target:Target, pattern:Option[IsotopePattern])
	
	def matchTargets(
			targets:Seq[Target], 
			patterns:Seq[IsotopePattern],
			specTime:Seq[Double]
	)(implicit params:DinosaurParams) = {
	  val apatterns = patterns.toArray
	  val sortedPatterns = apatterns.sortBy(_.hills.head.total.centerMz)
	  val patternMzs = sortedPatterns.map(_.hills.head.total.centerMz)
		for (t <- targets) yield {
		  if (params.verbose)
  			println("Matching target " + t.id)
  	  
  	  val (patternsStartIndx, patternsEndIndx) = DinoUtil.getMinxMaxIndx(patternMzs, t.mz, t.mzDiff)		  
			val candidates = for {
			  i <- patternsStartIndx until patternsEndIndx
			  val ip = sortedPatterns(i)
			  val rt = ip.apexHill.accurateApexRt(specTime)
			  if (rt > t.rtStart && rt < t.rtEnd && 
			        ip.apexHill.apex.intensity >= t.minApexInt && t.z == ip.z)
		  } yield ip
				
			val targetMidRt = (t.rtStart + t.rtEnd) / 2
			val best =
				if (candidates.nonEmpty) 
					Some(
						if (params.targetPreference.value == "close")
							candidates.minBy(p => math.abs(p.apexHill.accurateApexRt(specTime) - targetMidRt))
						else 
							candidates.maxBy(_.apexHill.apex.intensity)
					)
				else None
			Match(t, best)
		}
  }
}
