package se.lth.immun

import Deisotoper._
import DinoReport._

import se.lth.immun.mzml.ghost.GhostSpectrum

import java.io.File
import java.io.IOException
import java.awt.Color
import java.awt.Dimension
import java.awt.Color
import java.awt.Graphics2D
import java.awt.image.BufferedImage
import javax.imageio.ImageIO


object TargetReport {
	
	val mzMissingMarginLow = 1.0
	val mzMissingMarginHigh = 2.0
	val rtMargin = 0.5
					
	def createReports(
			targetMatches:Seq[TargetMatcher.Match],
			isotopes:Seq[IsotopePattern],
			hills:Seq[Hill],
			mzmlReader:MzMLReader, 
			streamer:ReportStreamer
	)(implicit params:DinosaurParams) = {
		
		val spectra = mzmlReader.specBackLog.reverse.toArray
		for (m <- targetMatches) {
			try {
				val targetIndex = mzmlReader.ms1IndexAtRt((m.target.rtEnd + m.target.rtStart)/2)
			
				m.pattern match {
					case Some(p) =>
						val targetReport = new TargetFoundReport(p, m.target, targetIndex, hills, spectra, params, mzmlReader.specTime, streamer)
						targetReport.report
						
					case None =>
						println("TARGETREPORT '%s': couldn't find any candidate isotope patterns".format(m.target.id))
						val relHills = relevantHills(m.target, hills, mzmlReader)
						val relIsotopes = isotopes.filter(iso => iso.hills.exists(h => relHills.contains(h)))
						val targetReport = new TargetMissingReport(
								m.target, targetIndex, 
								relHills,
								relIsotopes, spectra, params, streamer)
						targetReport.report
				}
			} catch {
				case e:Exception =>
					println("TARGETREPORT '%s': coultn't create report".format(m.target.id))
					if (params.verbose)
						e.printStackTrace()
			}
		}
	}
	
	
	def relevantHills(t:Target, hills:Seq[Hill], r:MzMLReader):Seq[Hill] = 
		hills.filter(h => {
			val mz = h.total.centerMz
			(mz >= t.mz - t.mzDiff - mzMissingMarginLow && mz < t.mz + t.mzDiff + mzMissingMarginHigh &&
			r.indexOverlapRt(h.scanIndex.head, h.scanIndex.last, t.rtStart - rtMargin, t.rtEnd + rtMargin))
		})
	
}


abstract class TargetReport(
		r:Target, 
		streamer:ReportStreamer,
		params:DinosaurParams
) extends DataViewable {
	
	/*def minInd:Int
	def maxInd:Int
	def minMz:Double
	def maxMz:Double
	
	def maxCPInt:Double
	def maxLogInt:Double
	
	val botPad = 3
	val hillH = 2
	def pw:Int
	def ph:Int
		
	def tox(ind:Int) = ((ind - minInd).toDouble / (maxInd - minInd) * pw).toInt
	def toy(mz:Double) = ph - ((mz - minMz) / (maxMz - minMz) * ph).toInt
	def toh(int:Double) = ((int / maxCPInt) * (50-botPad)).toInt
	def toLogCol(int:Double, minLogInt:Double, maxLogInt:Double):Color = {
		val relInt = math.max(0.0, (math.log10(1+int)-minLogInt)) / (maxLogInt-minLogInt)
		INT_GRADIENT.color(relInt)
	}
	def tocol(int:Double):Color = {
		val relInt = math.log10(1+int) / maxLogInt
		INT_GRADIENT.color(relInt)
	}
	
	def drawSpectra(
			g:Graphics2D, 
			relevantData:Seq[(Seq[(Double, Double)], Int)]
	):(Int, Int) = {
	
		val minInt = relevantData.flatMap(_._1.map(_._2)).filter(_>0.1).min
		val maxInt = relevantData.flatMap(_._1.map(_._2)).max
		val dw = math.min(10, math.max(1, pw / relevantData.length))
		val dh = math.min(10, math.max(1, ph / relevantData.map(_._1.length).max))
		
		val buffer = new Array[Array[Double]](pw)
		for (i <- 0 until pw) buffer(i) = new Array(ph)
		val bufferImage = new BufferedImage(pw, ph, BufferedImage.TYPE_INT_RGB)
		
		def updateBuffer(mz:Double, mzw:Double, index:Int, int:Double) = {
			val x1 = tox(index)
			val x2 = tox(index + 1)
			val y2 = toy(mz - mzw)
			val y1 = toy(mz + mzw)
			val w = (x2-x1)
			val h = (y2-y1)
			for {
				x <- x1 until math.min(pw, x1+math.max(w,1))
				y <- y1 until math.min(ph, y1+math.max(h,1))
			} {
				if (x < 0 || y < 0 || x >= pw || y >= ph)
					{}//println("ha")
				else
					buffer(x)(y) += int
			}
		}
		
		for ((data, i) <- relevantData) {
			val index = i + minInd
			if (data.nonEmpty) {
				for (((mz1, int),(mz2, _)) <- data.zip(data.tail))
					if (int > 0)
						updateBuffer(mz1, math.abs(mz2-mz1)/2, index, int)
				if (data.length == 1) {
					val (mz, int) = data.last
					if (int > 0)
						updateBuffer(mz, dh, index, int)
				}
				if (data.length > 1) {
					val (mz1, int) = data.last
					val (mz2, _) = data(data.length-2)
					if (int > 0)
						updateBuffer(mz1, math.abs(mz2-mz1)/2, index, int)
				}
			}
		}
		
		val bufferMin = math.log10(buffer.map(_.min).min+1)
		val bufferMax = math.log10(buffer.map(_.max).max+1)
		for {
			x <- 0 until pw
			y <- 0 until ph
		} 
			bufferImage.setRGB(x, y, toLogCol(buffer(x)(y), bufferMin, bufferMax).getRGB)
		
		g.drawImage(bufferImage, null, 0, 0)
		(dw, dh)
	}
	
	
	
	def drawPattern(g:Graphics2D, ip:IsotopePattern, col:Color, edgeText:Boolean) = {
		g.setColor(col)
		for (hill <- ip.hills)
			drawAnnotHill(g, hill)
		for ((h1, h2) <- ip.hills.zip(ip.hills.tail)) {
			g.setColor(col.darker())
			val x1 = tox(h1.scanIndex(h1.scanIndex.length / 2))
			val y1 = toy(h1.total.centerMz+h1.maxMzWidth/2)-1
			val x2 = tox(h2.scanIndex(h2.scanIndex.length / 2))
			val y2 = toy(h2.total.centerMz-h2.maxMzWidth/2)
			g.drawLine(x1, y1, x2, y2)
			
			if (edgeText) {
				val x = (x1+x2) / 2 + 5
				val y = (y1+y2) / 2 - 10
				val stats = new HillHillStats(h1, h2, ip.z, params)
				
				g.drawString("%.5f".format(stats.resMDiff), x, y)
				//g.drawString("%.5f".format(stats.totError), x, y + 15)
				g.drawString("x:%.3f".format(stats.corr), x, y + 15)
			}
			
			g.setColor(Color.WHITE)
			g.drawRect(x1-1, y1-1, 3, 3)
			g.drawRect(x2-1, y2-1, 3, 3)
		}
	}
	
	
	def drawSmallHill(g:Graphics2D, hill:Hill) = {
		val y1 = toy(hill.total.centerMz+hill.maxMzWidth/2)
		val y2 = toy(hill.total.centerMz-hill.maxMzWidth/2)
		g.drawRect(tox(hill.scanIndex.head), y1-1, 
				tox(hill.scanIndex.last+1)-tox(hill.scanIndex.head), (y2-y1)+1)
	}
	
	
	def drawAnnotHill(g:Graphics2D, hill:Hill) = {
		val y1 = toy(hill.total.centerMz+hill.maxMzWidth/2)
		val y2 = toy(hill.total.centerMz-hill.maxMzWidth/2)
		g.drawRect(tox(hill.scanIndex.head), y1-1, 
				tox(hill.scanIndex.last+1)-tox(hill.scanIndex.head), (y2-y1)+1)
		if (hill.scanIndex.length > 3)
			g.drawString("%.2f".format(hill.total.centerMz), tox(hill.scanIndex.head), y1-5)
	}
	
	def drawSplits(g:Graphics2D, hill:Hill) = 
		if (hill.leftSplit || hill.rightSplit) {
			val x1 = tox(hill.scanIndex.head)
			val x2 = tox(hill.scanIndex.last+1)
			val y1 = toy(hill.total.centerMz+hill.maxMzWidth/2)
			val y2 = toy(hill.total.centerMz-hill.maxMzWidth/2)
			if (hill.leftSplit) g.drawLine(x1, y1-1, x1, y2+1)
			if (hill.rightSplit) g.drawLine(x2, y1-1, x2, y2+1)
		}
	
	def drawHill(
			g:Graphics2D,
			hill:Hill
	) = {
		val mzw = hill.maxMzWidth //hill.maxMz.zip(hill.minMz).map(t => t._1 - t._2).max
		val minMz = hill.total.minMz - mzw/2//hill.minMz.min - mzw/2
		val maxMz = hill.total.maxMz + mzw/2 
		
		
	}
	*/
	def repType:String
	def render(g:Graphics2D):Unit
	def report = {
		val dir = r.group.map(_ + "/").getOrElse("")
		DinoReport.makeReport(
				streamer,
				"targets/%s%s%s.png".format(dir, r.id, repType),
				pw, ph,
				render _)
	}	
}

class TargetFoundReport(
		val ip:IsotopePattern,
		r:Target,
		targetIndex:Int,
		hills:Seq[Hill],
		allSpectra:Seq[GhostSpectrum],
		params:DinosaurParams,
		specTime:Seq[Double], 
		streamer:ReportStreamer
) extends TargetReport(r, streamer, params) {
	
	val minInd = ip.hills.map(_.scanIndex.min).min - 10
	val maxInd = ip.hills.map(_.scanIndex.max).max + 10
	val minMz = ip.hills.map(_.total.centerMz).min - 1
	val maxMz = ip.hills.map(_.total.centerMz).max + 1
	
	
	val relHills = hills.filter(h => {
			val mz = h.total.centerMz
			mz >= minMz && mz < maxMz &&
			h.scanIndex.last >= minInd && h.scanIndex.head < maxInd
		})
		
	
	val relSpectra = allSpectra.slice(minInd, maxInd)
	val relevantData = relSpectra.map(gs => 
		gs.mzs.zip(gs.intensities)
			.filter(t => t._1 > minMz && t._1 < maxMz)
			.sortBy(_._1)
		).zipWithIndex
	
	val maxInt = math.max(relevantData.flatMap(_._1.map(_._2)).max, params.adv.reportTargetIntMax)
	val maxLogInt = math.log10(1+maxInt)
	val maxCPInt = math.max(
			ip.hills.map(_.smoothIntensity.max).max, 
			ip.hills.map(_.rawIntensity.max).max)
	
	var pw = 1200
	var ph = 1000
	
	def repType:String = ".ok"
	def render(g:Graphics2D) = {
		val (dw, dh) = drawSpectra(g, relevantData)
		//g.setColor(Color.GRAY)
		//for (h <- relHills)
		//	drawHill(g, h)
		
		drawPattern(g, ip, Color.RED, pw/1000, true, true)(params)
		
		g.setColor(new Color(0.0f, 0.0f, 0.0f, 0.5f))
		g.fillRect(0, 0, pw, 50)
		
		g.setColor(Color.RED)
		g.drawString("Report target '%s':  mz=%.3f  rtStart=%.1f  rtEnd=%.1f  minApexInt=%.1e  apexIndex=%d".format(
				r.id, r.mz, r.rtStart, r.rtEnd, r.minApexInt, targetIndex
			), 10, 15)
		g.drawString("IsotopePattern  monoisoMz=%.3f  z=%d  apexInt=%.1e  apexRt=%.2f".format(
					ip.hills.head.total.centerMz, ip.z, ip.apexHill.apex.intensity, ip.apexHill.accurateApexRt(specTime)),
				10, 30)
		g.drawString("mz=[%.1f-%.1f]".format(minMz, maxMz), 10, 45)
		
		
		val tx = tox(targetIndex)
		val ty = toy(r.mz)
		g.setColor(Color.WHITE)
		g.drawOval(tx-15, ty-15, 30, 30)
	}
}


class TargetMissingReport(
		r:Target,
		targetIndex:Int,
		hills:Seq[Hill],
		isotopes:Seq[IsotopePattern],
		allSpectra:Seq[GhostSpectrum],
		params:DinosaurParams, 
		streamer:ReportStreamer
) extends TargetReport(r, streamer, params) {
	
	val minInd = hills.map(_.scanIndex.min).min
	val maxInd = hills.map(_.scanIndex.max).max
	val minMz = math.min(hills.map(_.total.centerMz).min, r.mz - r.mzDiff - 0.5)
	val maxMz = math.min(hills.map(_.total.centerMz).max, r.mz + r.mzDiff + 0.5)
	
	val relSpectra = allSpectra.slice(minInd, maxInd)
	val relevantData = relSpectra.map(gs => 
		gs.mzs.zip(gs.intensities)
			.filter(t => t._1 > minMz && t._1 < maxMz)
			.sortBy(_._1)
		).zipWithIndex
	
	val maxInt = relevantData.flatMap(_._1.map(_._2)).max
	val maxLogInt = math.log10(1+maxInt)
	val maxCPInt = math.max(
			hills.map(_.smoothIntensity.max).max, 
			hills.map(_.rawIntensity.max).max)
			
	var pw = 1200
	var ph = 1000
	
	def repType:String = ".miss"
	def render(g:Graphics2D) = {
		val (dw, dh) = drawSpectra(g, relevantData)
		
		g.setColor(Color.RED)
		for (h <- hills)
			drawHill(g, h, true, pw/1000)
		
		val sortedIsotopes = isotopes.sortBy(_.apexHill.apex.intensity).reverse
		for (iso <- sortedIsotopes.drop(10))
			drawPattern(g, iso, Color.DARK_GRAY, pw/1000, false, true)(params)
		for (iso <- sortedIsotopes.take(10))
			drawPattern(g, iso, Color.RED, pw/1000, true, true)(params)
		
			
		g.setColor(Color.WHITE)
		for (h <- hills)
			drawSplits(g, h)
		
		g.setColor(new Color(0.0f, 0.0f, 0.0f, 0.5f))
		g.fillRect(0, 0, pw, 35)
		g.setColor(Color.RED)
		g.drawString("Report target '%s':  mz=%.3f  rtStart=%.1f  rtEnd=%.1f  minApexInt=%.1e  apexIndex=%d".format(
				r.id, r.mz, r.rtStart, r.rtEnd, r.minApexInt, targetIndex
			), 10, 15)
		g.drawString("mz=[%.1f-%.1f]".format(minMz, maxMz), 10, 30)
		
		val tx = tox(targetIndex)
		val ty = toy(r.mz)
		g.setColor(Color.WHITE)
		g.drawOval(tx-15, ty-15, 30, 30)
		//g.drawLine(tx - 5, ty, tx + 5, ty)
		//g.drawLine(tx, ty - 5, tx, ty + 5)
	}	
}
	