package se.lth.immun


import java.io.File
import java.io.IOException
import java.awt.Color
import java.awt.Dimension
import java.awt.Graphics2D
import java.awt.image.BufferedImage
import javax.imageio.ImageIO

import scala.util.Try

import Deisotoper._

object NonLinMassCalibReport {

	val pw = 500
	val ph = 500
	

	def createReports(
			isotopes:Seq[IsotopePattern], 
			massCalib:Try[MzIntensity => Double], 
			streamer:ReportStreamer
	) = {
		try {
			NonLinMassCalibReport(isotopes, massCalib.get, streamer)
		} catch {
			case e:Exception =>
				println("MASSCALIBREPORT: failed to create report: " + e.getMessage)
		}
	}
	
	
	def apply(
			patterns:Seq[IsotopePattern], 
			f:MzIntensity => Double, 
			streamer:ReportStreamer
	) = {
		val massDiffs = patterns.map(p => (p.mass, (f(MzIntensity(p.mass, p.intensity)) - p.mass) / p.mass * 1e6))
		
		val minMass = massDiffs.map(_._1).min
		val maxMass = massDiffs.map(_._1).max
		val maxAbsDiff = math.max(massDiffs.map(_._2).max, math.abs(massDiffs.map(_._2).min))
		
		def tox(mass:Double) = ((mass - minMass) / (maxMass - minMass) * pw).toInt
		def toy(massDiff:Double) = ph - ((massDiff + maxAbsDiff) / (2*maxAbsDiff) * ph).toInt
	
		val nLines = 10
		
		def draw(g:Graphics2D) = {
			g.setColor(Color.WHITE)
			g.fillRect(0, 0, pw, ph)
			
			g.setColor(Color.LIGHT_GRAY)
			val fm = g.getFontMetrics

			for (mass <- (1 until nLines).map(i => minMass + (maxMass - minMass)/nLines*i)) {
				g.setColor(Color.LIGHT_GRAY)
				g.drawLine(tox(mass), 0, tox(mass), ph)
				val str = "%.1f".format(mass)
				val w = fm.charsWidth(str.toArray, 0, str.length)
				g.setColor(Color.GRAY)
				g.drawString(str, tox(mass) - w/2, ph - 5 - fm.getDescent)
			}
			for (diff <- (1 until nLines).map(i => -maxAbsDiff + 2*maxAbsDiff/nLines*i)) {
				g.setColor(Color.LIGHT_GRAY)
				g.drawLine(0, toy(diff), pw, toy(diff))
				val str = "%.2f".format(diff)
				val w = fm.charsWidth(str.toArray, 0, str.length)
				g.setColor(Color.GRAY)
				g.drawString(str, pw - w - 5, toy(diff) + fm.getDescent)
			}
			
			g.setColor(Color.BLUE)
			for ((mass, diff) <- massDiffs) 
				g.drawRect(tox(mass), toy(diff), 2, 2)
			
			g.setColor(Color.BLACK)
			g.drawString("mass correction, ppm vs mass, n="+patterns.length, 5, 5 + fm.getAscent)
				
		}
		
		DinoReport.makeReport(
				streamer,
				"nonLinMassCalibration.png",
				pw, ph,
				draw
			)
	}
}