package se.lth.immun

import java.io.File
import java.io.IOException
import java.awt.Color
import java.awt.Dimension
import java.awt.Graphics2D
import java.awt.image.BufferedImage
import javax.imageio.ImageIO

import se.lth.immun.mzml.ghost.GhostSpectrum

import scala.collection.mutable.HashMap
import DinoReport._

object HillReport {

	
	def createReports(
			hills:Seq[Hill], 
			specMap:HashMap[Int, Array[GhostSpectrum]], 
			streamer:ReportStreamer,
			params:DinosaurParams
	) = {
		for ((ms1start, spectra) <- specMap) {
			
			try {
				val reportHills = hills.filter(h => 
					h.scanIndex.head >= ms1start &&
					h.scanIndex.last < ms1start + spectra.length)
				if (reportHills.nonEmpty)
					HillReport(reportHills, ms1start, spectra, streamer, params)
			} catch {
				case e:Exception =>
					println("HILLREPORT %d-%d: failed to create report: ".format(ms1start, ms1start+spectra.length) + e.getMessage)
			}
		}
	}
	
	
	
	def apply(
			hills:Seq[Hill],
			ms1start:Int,
			spectra:Seq[GhostSpectrum],
			streamer:ReportStreamer,
			params:DinosaurParams
	):Unit = {
		if (spectra.isEmpty)
			throw new Exception("No spectra to build heatmap!")
	
		if (hills.isEmpty)
			throw new Exception("No hill for hill report!")
		
		val w = if (params.reportHighRes) 3000 else 400
		val h = if (params.reportHighRes) 4000 else 600
		val sorted = hills.sortBy(_.length).filter(_.length > 2)
				
		DinoReport.makeReport(
				streamer, 
				"hills_%d_%.3f.png".format(ms1start, sorted.last.total.centerMz),
				w, h,
				g => {
					
					if (params.reportHighRes) {
						g.setStroke(DinoReport.HIGHRES_STROKE)
						g.setFont(DinoReport.HIGHRES_FONT)
					}
					
					drawHill(sorted.last, ms1start, spectra, g, w, h/3, params)
					g.translate(0, h/3)
					drawHill(sorted(sorted.length / 2), ms1start, spectra, g, w, h/3, params)
					g.translate(0, h/3)
					drawHill(sorted.head, ms1start, spectra, g, w, h/3, params)
				})
	}
	
	
	
	
	def drawHill(
			hill:Hill,
			ms1start:Int,
			spectra:Seq[GhostSpectrum],
			g:Graphics2D,
			w:Int,
			h:Int,
			params:DinosaurParams
	) = {
		
		val minInd = hill.scanIndex.min - 10
		val maxInd = hill.scanIndex.max + 10
		val sMinInd = math.max(minInd, ms1start)
		val sMaxInd = math.min(maxInd, ms1start + spectra.length)
		
		val dMinMz = spectra.map(s => if (s.mzs.isEmpty) 0 else s.mzs.min).min
		val dMaxMz = spectra.map(s => if (s.mzs.isEmpty) 1 else s.mzs.max).max
		val mzw = hill.maxMzWidth 
		val minMz = math.max(dMinMz, hill.total.minMz - mzw/2)
		val maxMz = math.min(dMaxMz, hill.total.maxMz + mzw/2) 
		val heatHeight = h*3/4
		val barHeight = h - heatHeight
		
		val relevantSpectra = spectra.slice(sMinInd-ms1start, sMaxInd-ms1start)
		val relevantData = relevantSpectra.map(gs => gs.mzs.zip(gs.intensities).filter(t => t._1 >= minMz && t._1 <= maxMz)).zipWithIndex
		
		val maxInt = math.max(relevantData.flatMap(_._1.map(_._2)).max, params.adv.reportTargetIntMax)
		val maxLogInt = math.log10(1+maxInt)
		val maxCPInt = math.max(hill.smoothIntensity.max, hill.rawIntensity.max)
		
		val domain = DataView.domain(minInd, maxInd, minMz, maxMz, maxCPInt, maxLogInt) _
		val heatMap = domain(w, heatHeight)
		val barGraph = domain(w, barHeight)
		
		heatMap.drawSpectra(g, relevantData)
		
		
		g.setColor(Color.RED)
		heatMap.drawHillCenter(g, hill)
		
		g.setColor(Color.RED)
		g.drawString("mz=%.4f".format(hill.total.centerMz), 10, 15)
		g.drawString("maxInt=%e".format(maxInt), 10, 30)
		g.drawString(""+hill.scanIndex(hill.length / 2), 300, 15)
		
		g.translate(0, heatHeight)
		g.setColor(Color.GRAY)
		g.fillRect(0, 0, w, barHeight)
		barGraph.drawHillProfile(g, hill, Color.RED, Color.WHITE)
		g.translate(0, -heatHeight)
	}
}