package se.lth.immun

object TargetMatcher {
	case class Match(target:Target, pattern:Option[IsotopePattern])
	
	def matchTargets(
			targets:Seq[Target], 
			patterns:Seq[IsotopePattern],
			specTime:Seq[Double]
	)(implicit params:DinosaurParams) = 
		for (t <- targets) yield {
			val candidates = 
				patterns.filter(ip => {
					val rt = ip.apexHill.accurateApexRt(specTime)
					val mz = ip.hills.head.total.centerMz
					rt >= t.rtStart && rt < t.rtEnd && 
					mz >= t.mz - t.mzDiff && mz < t.mz + t.mzDiff &&
					ip.apexHill.apex.intensity >= t.minApexInt &&
					t.z == ip.z
				})
				
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