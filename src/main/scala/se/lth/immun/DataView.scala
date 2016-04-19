package se.lth.immun

import java.awt.Color
import java.awt.Graphics2D
import java.awt.image.BufferedImage
import DinoReport._


object DataView {
	def domain(
		minInd:Int,
		maxInd:Int,
		minMz:Double,
		maxMz:Double,
		maxCPInt:Double,
		maxLogInt:Double
	)(
		pw:Int,
		ph:Int
	) = new DataView(minInd, maxInd, minMz, maxMz, maxCPInt, maxLogInt, pw, ph)
}

class DataView(
		val minInd:Int,
		val maxInd:Int,
		val minMz:Double,
		val maxMz:Double,
		val maxCPInt:Double,
		val maxLogInt:Double,
		val pw:Int,
		val ph:Int
) extends DataViewable 


trait DataViewable {

	def minInd:Int
	def maxInd:Int
	def minMz:Double
	def maxMz:Double
	
	def maxCPInt:Double
	
	val botPad = 3
	val hillH = 2
	def pw:Int
	def ph:Int
		
	def tox(ind:Int) = ((ind - minInd).toDouble / (maxInd - minInd) * pw).toInt
	def toy(mz:Double) = ph - ((mz - minMz) / (maxMz - minMz) * ph).toInt
	def toh(int:Double) = ((int / maxCPInt) * (ph-botPad)).toInt
	def toLogCol(int:Double, minLogInt:Double, maxLogInt:Double):Color = {
		val relInt = math.max(0.0, (math.log10(1+int)-minLogInt)) / (maxLogInt-minLogInt)
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
	
	
	
	def drawPattern(
			g:Graphics2D, 
			ip:IsotopePattern, 
			col:Color, 
			lineWidth:Int,
			edgeText:Boolean,
			hillText:Boolean
	)(
			implicit params:DinosaurParams
	) = {
		g.setColor(col)
		for (hill <- ip.hills)
			drawHill(g, hill, hillText, lineWidth)
		for ((h1, h2) <- ip.hills.zip(ip.hills.tail)) {
			g.setColor(col.darker())
			val x1 = tox(h1.scanIndex(h1.scanIndex.length / 2))
			val y1 = toy(h1.total.centerMz+h1.maxMzWidth/2)-lineWidth
			val x2 = tox(h2.scanIndex(h2.scanIndex.length / 2))
			val y2 = toy(h2.total.centerMz-h2.maxMzWidth/2)
			g.drawLine(x1, y1, x2, y2)
			
			if (edgeText) {
				val x = (x1+x2) / 2 + 5*lineWidth
				val y = (y1+y2) / 2 - 10*lineWidth
				val stats = new HillHillStats(h1, h2, ip.z, params)
				
				g.drawString("%.5f".format(stats.resMDiff), x, y)
				g.drawString("x:%.3f".format(stats.corr), x, y + 15*lineWidth)
			}
			
			g.setColor(Color.WHITE)
			g.drawRect(x1-1, y1-1, 3, 3)
			g.drawRect(x2-1, y2-1, 3, 3)
		}
	}
	
	
	def drawHill(g:Graphics2D, hill:Hill, annot:Boolean, lineWidth:Int) = {
		val y1 = toy(hill.total.centerMz+hill.maxMzWidth/2)
		val y2 = toy(hill.total.centerMz-hill.maxMzWidth/2)
		g.drawRect(tox(hill.scanIndex.head), y1-lineWidth, 
				tox(hill.scanIndex.last+lineWidth)-tox(hill.scanIndex.head), (y2-y1)+lineWidth)
		if (annot && hill.scanIndex.length > 3)
			g.drawString("%.2f".format(hill.total.centerMz), tox(hill.scanIndex.head), y1-5*lineWidth)
	}
	
	
	def drawHillCenter(g:Graphics2D, hill:Hill) = {
		val bins = 
			for (i <- 0 until hill.length) yield {
				val scan = hill.scanIndex(i)
				val x0 = tox(scan)
					//if (scan <= minInd) tox(minInd)
					//else (tox(scan-1) + tox(scan)) / 2
				val x1 = tox(scan+1)
					//if (scan >= maxInd) tox(maxInd)
					//else (tox(scan) + tox(scan+1)) / 2
				(x0, x1)
			}
		for (i <- 0 until hill.length) {
			val (x0, x1) = bins(i)
			val y = toy(math.min(maxMz, math.max(minMz, hill.centerMz(i))))
			g.drawLine(x0, y, x1, y)
		}
		for (i <- 1 until hill.length) {
			val (_, x0) = bins(i-1)
			val (x1, _) = bins(i)
			val y0 = toy(math.min(maxMz, math.max(minMz, hill.centerMz(i-1))))
			val y1 = toy(math.min(maxMz, math.max(minMz, hill.centerMz(i))))
			g.drawLine(x0, y0, x1, y1)
		}
	}
	
	
	def drawHillProfile(g:Graphics2D, hill:Hill, colRaw:Color, colSmooth:Color) = {
		val bins = 
			for (i <- 0 until hill.length) yield {
				val scan = hill.scanIndex(i)
				val x0 = 
					if (scan <= minInd) tox(minInd)
					else (tox(scan-1) + tox(scan)) / 2
				val x1 = 
					if (scan >= maxInd) tox(maxInd)
					else (tox(scan) + tox(scan+1)) / 2
				(x0, x1)
			}
		g.setColor(colRaw)
		for (i <- 0 until hill.length) {
			val (x0, x1) = bins(i)
			val h = toh(hill.rawIntensity(i))
			g.fillRect(x0, ph-botPad-h, x1-x0, h)
		}
		g.setColor(colSmooth)
		for (i <- 1 until hill.length) {
			val x0 = tox(hill.scanIndex(i-1))
			val x1 = tox(hill.scanIndex(i))
			val h0 = toh(hill.smoothIntensity(i-1))
			val h1 = toh(hill.smoothIntensity(i))
			g.drawLine(x0, ph-botPad-h0, x1, ph-botPad-h1)
		}
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
}