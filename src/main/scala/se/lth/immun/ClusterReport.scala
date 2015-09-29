package se.lth.immun

import java.io.File
import java.io.IOException
import java.awt.Color
import java.awt.Dimension
import java.awt.Graphics2D
import java.awt.image.BufferedImage
import javax.imageio.ImageIO

import se.lth.immun.mzml.ghost.GhostSpectrum
import DinoReport._
import Cluster.Edge

class ClusterReport(
		val clusters:Seq[Seq[Edge]],
		val hills:Seq[Hill],
		val spectra:Seq[GhostSpectrum],
		val params:DinosaurParams
) {
	
	implicit val dp = params
	val specMinInd = spectra.map(_.spectrum.index).min
	val specMaxInd = spectra.map(_.spectrum.index).max
	
	val cs = clusters.map(c => (Cluster(c, hills), c))
	
	val inSpectra = cs.filter(t => 
		t._1.hills.map(_.scanIndex.min).min < specMaxInd && 
		t._1.hills.map(_.scanIndex.max).max > specMinInd)
	
	val maxMzHeight = 10
	val minMz = math.floor(inSpectra.map(_._1.hills.map(_.total.centerMz).min).min)
	val maxMz = math.ceil(math.min(minMz + maxMzHeight, inSpectra.map(_._1.hills.map(_.total.centerMz).max).max))
	
	val visible = inSpectra.filter(t => t._1.hills.map(_.total.centerMz).max < maxMz)
	
	val pw = 1000
	val ph = 1000
	
	def tox(ind:Int) = ((ind - specMinInd).toDouble / (specMaxInd - specMinInd) * pw).toInt
	def toy(mz:Double) = ph - ((mz - minMz) / (maxMz - minMz) * ph).toInt
	def toLogCol(int:Double, minLogInt:Double, maxLogInt:Double):Color = {
		val relInt = math.max(0.0, (math.log10(1+int)-minLogInt)) / (maxLogInt-minLogInt)
		INT_GRADIENT.color(relInt)
	}
	
	
	def report = {
		val image = new BufferedImage(1000, 1000, BufferedImage.TYPE_INT_RGB)
		val g = image.createGraphics
		
		val (dw, dh) = drawSpectra(g, spectra)
		for ((c, es) <- visible)
			drawCluster(g, c, dw)
		
		g.setColor(Color.RED)
		g.drawString("mz=[%.1f-%.1f]".format(minMz, maxMz), 10, 15)
		g.drawString("scans=[%d-%d]".format(specMinInd, specMaxInd), 10, 30)
		
		g.dispose
		try { 
		    ImageIO.write(image, "png", new File("qc/clusters.png")) 
		} catch {
			case ioe:IOException =>
		    	ioe.printStackTrace
		}
	}
	
	
	
	def drawSpectra(g:Graphics2D, spectra:Seq[GhostSpectrum]):(Int, Int) = {
		
		val relevantData = spectra.map(gs => 
				(gs.spectrum.index, 
				 gs.mzs.zip(gs.intensities).filter(t => t._1 > minMz && t._1 < maxMz)
			))
		
		val minInt = relevantData.map(_._2.map(_._2).filter(_>0.1).min).min
		val minLogInt = math.log10(1+minInt)
		val maxInt = relevantData.map(_._2.map(_._2).max).max
		val maxLogInt = math.log10(1+maxInt)
		val dw = math.min(10, math.max(1, pw / relevantData.length))
		val dh = math.min(10, math.max(1, ph / relevantData.map(_._2.length).max))
		
		for {	
			(index, data) <- relevantData
			(mz, int) <- data
		} {
			if (int > 0) {
				val x = tox(index)
				val y = toy(mz)
				g.setColor(toLogCol(int, minLogInt, maxLogInt))
				g.fillRect(x-dw, y-dh, dw, dh)
			}
		}
		
		(dw, dh)
	}
	
	val hillH = 2
	var i = 0
	def drawCluster(g:Graphics2D, c:Cluster, dw:Int) = {
		g.setColor(RED4.color(i))
		i += 1
		for (hill <- c.hills)
			drawHill(g, hill, dw)
	}
	
	
	
	def drawHill(g:Graphics2D, hill:Hill, dw:Int) = {
		
		g.drawRect(tox(hill.scanIndex.head), toy(hill.total.centerMz)-hillH-1, 
				tox(hill.scanIndex.last)-tox(hill.scanIndex.head), 2*hillH+1)
		/*for ((ind, mz) <- hill.scanIndex.zip(hill.centerMz)) 
			g.drawLine(tox(ind)-dw, toy(mz), tox(ind)+dw, toy(mz))
			
		val indMz = hill.scanIndex.zip(hill.centerMz)
		for (((pInd,pMz), (nInd,nMz)) <- indMz.zip(indMz.tail))
			g.drawLine(tox(pInd)-dw, toy(pMz), tox(nInd)+dw, toy(nMz))
			* */
	}

}