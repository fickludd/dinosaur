package se.lth.immun

import java.io.File
import java.io.IOException
import java.awt.Color
import java.awt.Dimension
import java.awt.Graphics2D
import java.awt.image.BufferedImage
import javax.imageio.ImageIO

import se.lth.immun.graphs.LineGraph
import se.lth.immun.graphs.util._
import se.lth.immun.chem._

object IsotopeDistributionReport {

	val w = 400
	val h = 200
	
	
	def createReports(
			isotopes:Seq[IsotopePattern],
			streamer:ReportStreamer
	)(implicit params:DinosaurParams) = {
		val shuffled1000 = params.reportRandom.shuffle(isotopes.zipWithIndex.take(1000))
		for ((ip, i) <- shuffled1000.take(params.nReport))
			try {
				IsotopeDistributionReport(ip, i, streamer)
			} catch {
				case e:Exception =>
					println("ISODISTREPORT %d: failed to create report: ".format(i) + e.getMessage)
			}
	}
	
	
	def apply(
			ip:IsotopePattern, 
			i:Int,
			streamer:ReportStreamer
	) = {
		DinoReport.makeReport(
				streamer,
				"ipattern_"+i+".png",
				w, h,
				g => {
					val ipGraph = patternGraph(g, ip)
					ipGraph.renderer.setup(
							ipGraph.xAxis, 
							ipGraph.yAxis, 
							ipGraph.style, 
							new Size(w, h)
						)
					ipGraph.render(g, ipGraph.renderer)
				})
	}
	
	def patternGraph(g:Graphics2D, ip:IsotopePattern) = {
		val profile = ip.hills.map(_.apex.intensity)
		val distr = Peptide.averagine(ip.mass).getIsotopeDistribution()
			
		val fg 	= new LineGraph
		fg.preferredSize 			= new Dimension(w, h)
		fg.style.annotColBackground = new Color(0.5f, 0.5f, 0.5f, 0.5f)
		fg.setCurves(List(
			new Curve2(
				ip.hills.map(_.total.centerMz),
				profile.map(_ / profile.max),
				profile.map(_ => false),
				"measured",
				Color.BLUE
			),
			new Curve2(
				distr.intensities.indices.map(i => 
					(ip.mass + (i+ip.isotopeOffset)*DinoUtil.ISOTOPE_PATTERN_DIFF) / 
						ip.z.toDouble + Constants.PROTON_WEIGHT
				),
				distr.intensities,
				distr.intensities.map(_ => false),
				"averagine",
				Color.GREEN
			)))
		fg.title = "offset=%d mass=%.1f corr=%.2f".format(
				ip.isotopeOffset, ip.mass, ip.averagineCorr
			)
		fg.repaint
		fg
	}
}