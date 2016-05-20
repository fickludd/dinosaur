package se.lth.immun

import java.io.Writer
import java.io.File
import Deisotoper._
import se.lth.immun.chem.Constants

import se.lth.immun.protocol.MsFeatureFile
import se.lth.immun.protocol.MsFeatures
import se.lth.immun.protocol.MsFeature
import se.lth.immun.protocol.MsHill
import se.lth.immun.protocol.MSFeatureProtocol.RtUnit

object OutWriter {
	def writeBinary(
			f:File, 
			features:Seq[IsotopePattern], 
			rtMap:Seq[Double],
			params:DinosaurParams
	):Long = {
		val t0 = System.currentTimeMillis
		MsFeatureFile.write(
			f, 
			MsFeatures(
				rtMap, 
				features.zipWithIndex.map(t => {
					val (p, i) = t
					MsFeature(
							i,
							p.hills.head.total.centerMz,
							p.z,
							p.mass,
							p.apexHill.accurateApexRt(rtMap),
							p.startRt(rtMap),
							p.endRt(rtMap),
							p.apexHill.apex.intensity,
							p.intensity,
							Some(p.averagineCorr),
							p.hills.map(h =>
								MsHill(
									h.scanIndex.head,
									h.scanIndex.last,
									h.total.centerMz,
									h.total.centerMzError,
									Some(h.fwhm(rtMap)),
									h.accurateApexRt(rtMap),
									h.apex.intensity,
									h.smoothFullProfile
							))
					)}
				)),
			RtUnit.MINUTE,
			false) // params.verbose) <- WAY TOO MUCH OUTPUT
		System.currentTimeMillis() - t0
	}
}

class OutWriter(val params:DinosaurParams, w:Writer) {
	
	def writeRow(qoute:Boolean)(a:Any*) = 
		w.write(a.map(_ match {
			case s:String => 
				if (qoute) params.outQuote + s + params.outQuote
				else s
			case x => x.toString
		}).mkString(params.outSep) + "\n")
	
		
		
	def writeFeatureCsv(
			result:DinosaurResult, 
			specTime:Seq[Double]
	):Long = {
		val t = System.currentTimeMillis
		
		writeRow(false)("mz", "mostAbundantMz", "charge", "rtStart", "rtApex", "rtEnd", "fwhm", "nIsotopes", "nScans", "averagineCorr", "mass", "massCalib", "intensityApex", "intensitySum")
		
		for (p <- result.patterns)
			writeRow(true)(
					p.hills.head.total.centerMz,
					p.hills(p.mostAbundNbr).total.centerMz,
					p.z,
					p.startRt(specTime),
					p.apexHill.accurateApexRt(specTime),
					p.endRt(specTime),
					p.apexHill.fwhm(specTime),
					p.hills.length,
					p.hills.map(_.scanIndex.length).max,
					p.averagineCorr,
					p.mass,
					result.massCalib(MzIntensity(p.mass, p.intensity)),
					p.apexHill.apex.intensity,
					p.intensity
			)
		
		w.close
		
		System.currentTimeMillis - t
	}
	
		
		
	def writeTargetCsv(
			result:DinosaurResult, 
			specTime:Seq[Double]
	):Long = {
		val t = System.currentTimeMillis
		
		writeRow(false)("mz", "mostAbundantMz", "charge", "rtStart", "rtApex", "rtEnd", "fwhm", "nIsotopes", "nScans", "averagineCorr", "mass", "massCalib", "intensityApex", "intensitySum", "targetGroup", "targetId")
		
		for (TargetMatcher.Match(t, pOpt) <- result.targetMatches)
			pOpt match {
				case None =>
					writeRow(true)(
						0.0,
						0.0,
						t.z,
						0.0,
						0.0,
						0.0,
						0.0,
						0,
						0,
						0.0,
						0.0,
						0.0,
						0,
						0,
						t.group,
						t.id
					)
				case Some(p) =>
					writeRow(true)(
						p.hills.head.total.centerMz,
						p.hills(p.mostAbundNbr).total.centerMz,
						p.z,
						p.startRt(specTime),
						p.apexHill.accurateApexRt(specTime),
						p.endRt(specTime),
						p.apexHill.fwhm(specTime),
						p.hills.length,
						p.hills.map(_.scanIndex.length).max,
						p.averagineCorr,
						p.mass,
						result.massCalib(MzIntensity(p.mass, p.intensity)),
						p.apexHill.apex.intensity,
						p.intensity,
						t.group,
						t.id
					)
			}
			
		
		w.close
		
		System.currentTimeMillis - t
	}
	
		
		
	def writeHillCsv(
			result:DinosaurResult,
			specTime:Seq[Double]
	):Long = {
		val t = System.currentTimeMillis
		
		writeRow(false)("patternId", "monoisoMz", "charge", "rtStart", "rtApex", "rtEnd", "fwhm", "nScans", "intensityApex", "intensitySum")
		
		for {
			(p, i) <- result.patterns.zipWithIndex
			h <- p.hills
		}
			writeRow(true)(
					i,
					h.total.centerMz,
					p.z, 
					h.startRt(specTime),
					h.accurateApexRt(specTime),
					h.endRt(specTime),
					h.fwhm(specTime),
					h.scanIndex.length,
					h.apex.intensity,
					h.total.intensity
			)
		
		w.close
		
		System.currentTimeMillis - t
	}
	
	
	/**
    	File adapter for MsInspect files.

		Lines with "#" are comments and are ignored.

		The first non-comment line is the header and contains the column names:<br>
		scan	time	mz	accurateMZ	mass	intensity	charge	chargeStates	kl	background	median	peaks	scanFirst	scanLast	scanCount	totalIntensity	sumSquaresDist	description

		Every subsequent line is a feature.
		
		OpenMS only uses cols 1,2,5,6 and 8. The rest are considered meta-data.
	 */
	def writeMsInspectCsv(
			result:DinosaurResult,
			specTime:Seq[Double]
	):Long = {
		val t = System.currentTimeMillis
		
		w.write("# file format intended for conversion to featureXML through OpenMS FileConverter\n")
		
		w.write("scan	time	mz	accurateMZ	mass	intensity	charge	chargeStates	kl	background	median	peaks	scanFirst	scanLast	scanCount	totalIntensity	sumSquaresDist	description\n")
		
		for {
			(p, i) <- result.patterns.zipWithIndex
		}
			writeRow(true)(
					i,
					p.apexHill.accurateApexRt(specTime) * 60,
					p.hills.head.total.centerMz,
					(p.mass + Constants.PROTON_WEIGHT * p.z) / p.z,
					p.mass,
					p.intensity,
					p.z,
					0,
					1.0,
					0.0,
					0.0,
					p.hills.length,
					p.startInd,
					p.endInd,
					p.hills.map(_.scanIndex.length).max,
					p.intensity,
					0.0,
					""
			)
			
		w.close
		
		System.currentTimeMillis - t
		
	}
	

}