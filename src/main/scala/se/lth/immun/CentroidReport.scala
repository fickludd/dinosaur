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

object CentroidReport {

	import Dinosaur._
	
	def apply(
			specIndex:Int, 
			mzs:Seq[Double], 
			ints:Seq[Double], 
			peaks:Seq[CentroidPeak], 
			streamer:ReportStreamer
	) = {
		val sorted = peaks.sortBy(_.int)
		val min = peakGraph(sorted.head, mzs, ints)
		val mid = peakGraph(sorted(sorted.length / 2), mzs, ints)
		val max = peakGraph(sorted.last, mzs, ints)
		
		val w = 600
		val h = 200
		
		DinoReport.makeReport(
				streamer,
				"centroid_"+specIndex+".png",
				w, h,
				g => {
					render(g, min, 0)
					render(g, mid, 200)
					render(g, max, 400)
				})
	}
	
	
	def render(g:Graphics2D, lg:LineGraph, x:Int) = {
		g.translate(x, 0)
		lg.renderer.setup(lg.xAxis, lg.yAxis, lg.style, new Size(200, 200))
		lg.render(g, lg.renderer)
		g.translate(-x, 0)
	}
	
	
	def peakGraph(p:CentroidPeak, mzs:Seq[Double], ints:Seq[Double]) = {
		val mzw = p.maxMz - p.minMz
		val minInd = mzs.indexWhere(_ > p.minMz - mzw)
		val maxInd = mzs.indexWhere(_ > p.maxMz + mzw)
				
		val fg 	= new LineGraph
		fg.preferredSize 			= new Dimension(200, 200)
		fg.style.annotColBackground = new Color(0.5f, 0.5f, 0.5f, 0.5f)
		fg.setCurves(List(
			new Curve2(
				mzs.slice(minInd, maxInd),
				ints.slice(minInd, maxInd),
				mzs.slice(minInd, maxInd).map(_ => false)
			)))
		fg.title = p.mz + "mz"
		fg.addAnnotation(new XAnnotation(p.minMz, Color.RED, "minMz"))
		fg.addAnnotation(new XAnnotation(p.mz, Color.RED, "mz"))
		fg.addAnnotation(new XAnnotation(p.maxMz, Color.RED, "maxMz"))
		fg.repaint
		fg
	}
}