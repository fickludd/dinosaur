package se.lth.immun

import java.io.File
import java.io.IOException
import java.awt.Color
import java.awt.Dimension
import java.awt.Graphics2D
import java.awt.image.BufferedImage
import javax.imageio.ImageIO

import scala.collection.mutable.HashMap
import se.lth.immun.mzml.ghost.GhostSpectrum
import Deisotoper._
import DinoReport._


object DeisotopeReport {
	
	def createReports(
			isotopes:Seq[IsotopePattern], 
			reader:MzMLReader,
			streamer:ReportStreamer
	)(implicit params:DinosaurParams) = {
		for ((ms1start, spectra) <- reader.ghostSpecMap) {
			try {
				val reportPatterns = isotopes.filter(ip => 
					ip.hills.map(_.scanIndex.min).min < ms1start + spectra.length && 
					ip.hills.map(_.scanIndex.max).max > ms1start)
				if (reportPatterns.isEmpty) { 
					println("DEISOREPORT %d-%d: no isotope patterns in spectrum range ".format(ms1start, ms1start+spectra.length))
				} else {
					val deisoReport = new DeisotopeReport(reportPatterns, ms1start, spectra, params)
					deisoReport.report(streamer)
				}
			} catch {
				case e:Exception =>
					println("DEISOREPORT %d-%d: failed to create report: ".format(ms1start, ms1start+spectra.length) + e.getMessage)
			}
		}
	}
}


class DeisotopeReport(
		val ipatterns:Seq[IsotopePattern],
		val ms1start:Int,
		val spectra:Seq[GhostSpectrum],
		params:DinosaurParams
) extends DataViewable {
	if (spectra.isEmpty)
		throw new Exception("No spectra to build heatmap!")

	val minInd = ms1start
	val maxInd = ms1start + spectra.length
	
	val mzHeight = params.reportDeisoMzHeight
	val dMinMz = spectra.map(s => if (s.mzs.isEmpty) 0 else s.mzs.min).min
	val dMaxMz = spectra.map(s => if (s.mzs.isEmpty) 1 else s.mzs.max).max
	val centerPattern = ipatterns(params.reportRandom.nextInt(ipatterns.length))
	val minMz = math.max(dMinMz, centerPattern.apexHill.total.centerMz - params.reportDeisoMzHeight * 0.33)
	val maxMz = math.min(dMaxMz, centerPattern.apexHill.total.centerMz + params.reportDeisoMzHeight * 0.67)
	
	val maxCPInt = 0.0
	
	val visible = ipatterns.filter(ip => 
			ip.hills.map(_.total.centerMz).min > minMz && 
			ip.hills.map(_.total.centerMz).max < maxMz)
	
	val pw = if (params.reportHighRes) 4000 else 1000
	val ph = if (params.reportHighRes) 4000 else 1000
	
	
	def report(streamer:ReportStreamer) = {
		
		var i = 0
		def nextColor = {
			i = (i + 1) % DinoReport.RED4.curveColors.length
			DinoReport.RED4.color(i)
		}
		
		DinoReport.makeReport(
				streamer,
				"deisotoper-%d.png".format(minInd),
				pw, ph,
				g => {
					
					if (params.reportHighRes) {
						g.setStroke(DinoReport.HIGHRES_STROKE)
						g.setFont(DinoReport.HIGHRES_FONT)
					}
					
					val relevantData = spectra.map(gs => 
							gs.mzs.zip(gs.intensities)
								.filter(t => t._1 > minMz && t._1 < maxMz)
						).zipWithIndex
					val (dw, dh) = drawSpectra(g, relevantData)
					
					val sorted = visible.sortBy(_.apexHill.apex.intensity).reverse
					for (ip <- sorted.drop(20)) 
						drawPattern(g, ip, nextColor, pw/1000, false, false)(params)
						
					for (ip <- sorted.take(20)) 
						drawPattern(g, ip, nextColor, pw/1000, false, true)(params)
					
					g.setColor(Color.RED)
					g.drawString("mz=[%.1f-%.1f]".format(minMz, maxMz), pw/100, ph*15/1000)
					g.drawString("scans=[%d-%d]".format(minInd, maxInd), pw/100, ph*30/1000)
				})
		
	}
}