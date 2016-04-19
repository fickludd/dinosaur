package se.lth.immun

import java.awt.Color
import java.awt.Graphics2D
import java.awt.image.BufferedImage
import java.awt.BasicStroke
import java.awt.Font
import javax.imageio.ImageIO
import java.io.IOException

import scala.util.Random

object DinoReport {

	class ColorGradient(
			val anchors:Seq[Double],
			val colors:Seq[Color]
	) {
		def color(x:Double):Color = {
			for (i <- 1 until anchors.length) {
				if (x >= anchors(i-1) && x <= anchors(i)) {
					val k = (x - anchors(i-1)) / (anchors(i) - anchors(i-1))
					return new Color(
							(colors(i-1).getRed * (1-k) + colors(i).getRed * k).toInt,
							(colors(i-1).getGreen * (1-k) + colors(i).getGreen * k).toInt,
							(colors(i-1).getBlue * (1-k) + colors(i).getBlue * k).toInt)
				}
			}
			throw new Exception(x+" outside gradient domain [%.2f,%.2f]".format(anchors.head, anchors.last))
		}
	}
	
	val INT_GRADIENT = new ColorGradient(Array(0.0, 0.35, 0.6, 0.75, 0.9, 1.0), Array(Color.BLACK, new Color(85, 0, 255), new Color(0, 85, 255), Color.YELLOW, new Color(0, 255, 230), Color.WHITE))
	val INT_GRADIENT2 = new ColorGradient(Array(0.0, 0.2, 0.6, 1.0), Array(Color.BLACK, Color.BLACK, Color.BLUE, Color.CYAN))
	val MONOCHROME_GRADIENT = new ColorGradient(Array(0.0, 1.0), Array(Color.WHITE, Color.BLACK))
	val HIGHRES_LINEWIDTH = 4.0
	val HIGHRES_STROKE = new BasicStroke((HIGHRES_LINEWIDTH-1).toFloat)
	val HIGHRES_FONT = new Font("TimesRoman", Font.PLAIN, 28)
	
	
	
	def randomColor = 
		new Color(
				(0.4 + Random.nextDouble*0.8).toFloat, 
				(0.4 + Random.nextDouble*0.8).toFloat, 
				(0.4 + Random.nextDouble*0.8).toFloat)
	
	def randomColorChannel(r:Boolean, g:Boolean, b:Boolean) = 
		new Color(
				if (r) (0.4 + Random.nextDouble*0.4).toFloat else 0.0f, 
				if (g) (0.4 + Random.nextDouble*0.4).toFloat else 0.0f, 
				if (b) (0.4 + Random.nextDouble*0.4).toFloat else 0.0f
			)
	
	
	/*
	 * Color taken from d3js category20
	 */
	val CATEGORY20 = new DistinctColors(Array(
			new Color(0x1f77b4), new Color(0xaec7e8), 
			new Color(0xff7f0e), new Color(0xffbb78), 
			new Color(0x2ca02c), new Color(0x98df8a), 
			new Color(0xd62728), new Color(0xff9896), 
			new Color(0x9467bd), new Color(0xc5b0d5), 
			new Color(0x8c564b), new Color(0xc49c94), 
			new Color(0xe377c2), new Color(0xf7b6d2), 
			new Color(0x7f7f7f), new Color(0xc7c7c7), 
			new Color(0xbcbd22), new Color(0xdbdb8d), 
			new Color(0x17becf), new Color(0x9edae5)))
	
	val RED4 = new DistinctColors(Array(Color.WHITE.darker, Color.RED, Color.MAGENTA, new Color(1.0f, 0.8f, 0.2f)))
	
	
	class DistinctColors(val curveColors:Seq[Color]) {
		val curveColorOrder:Seq[Int] = Random.shuffle((0 until curveColors.length).toArray)
		def color(i:Int) = curveColors(curveColorOrder(i % curveColors.length))
	}
	
	
	
	def makeReport(
			streamer:ReportStreamer, 
			relPath:String, 
			w:Int, h:Int, 
			f:Graphics2D => Unit,
			pdf:Boolean = false
	) = {
		
		if (pdf) {
			
		} else {
			val image = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB)
			val g = image.createGraphics
			
			f(g)
			
			g.dispose
			try { 
			    ImageIO.write(image, "png", streamer.streamByPath(relPath)) 
			} catch {
				case ioe:IOException =>
			    	ioe.printStackTrace
			}
		}
		streamer.closeLast
	}
}