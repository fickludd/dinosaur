package se.lth.immun

class Timer {
	var t0 = 0L
	def start =
		t0 = System.currentTimeMillis
		
	def click = {
		val now = System.currentTimeMillis
		val dt = now - t0
		t0 = now
		dt
	}
	
	def check = System.currentTimeMillis - t0
}

trait Timeable {

	val timer = new Timer
}