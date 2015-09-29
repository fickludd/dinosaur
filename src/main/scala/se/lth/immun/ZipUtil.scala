package se.lth.immun

import java.io.File
import java.io.IOException
import java.io.InputStream
import java.io.OutputStream
import java.io.FileOutputStream
import java.io.FileInputStream
import java.util.zip.ZipOutputStream
import java.util.zip.ZipEntry

object ZipUtil {

    def zipFolder(folder:File, zipFile:File, verbose:Boolean = false):Unit = 
        zipFolder(folder, new FileOutputStream(zipFile), zipFile.getName.dropRight(4), verbose)
    
    def zipFolder(folder:File, outputStream:OutputStream, outBaseName:String, verbose:Boolean):Unit = {
        val zipOutputStream = new ZipOutputStream(outputStream)
        processFolder(folder.getAbsoluteFile, zipOutputStream, outBaseName, folder.getAbsolutePath.length() + 1, verbose)
        zipOutputStream.close
    }

    def processFolder(
    		folder:File, 
    		zipOutputStream:ZipOutputStream, 
    		outBaseName:String, 
    		prefixLength:Int, 
    		verbose:Boolean
    ):Unit = 
        for (file <- folder.listFiles) {
            if (file.isFile()) {
                val zipEntry = new ZipEntry(outBaseName+"/"+file.getPath.substring(prefixLength))
                if (verbose)
                	println("adding entry '%s' to zip archive".format(zipEntry.getName))
                zipOutputStream.putNextEntry(zipEntry)
                val inputStream = new FileInputStream(file)
                streamCopy(inputStream, zipOutputStream)
                inputStream.close
                zipOutputStream.closeEntry()
            } else if (file.isDirectory)
                processFolder(file, zipOutputStream, outBaseName, prefixLength, verbose)
            
        }
    
    def streamCopy(input:InputStream, output:OutputStream) = {
    	val buffer = new Array[Byte](1024 * 4)
        var read = input.read(buffer)
        while (read != -1) { 
        	 output.write(buffer, 0, read); 
        	 read = input.read(buffer)
        }
    }
}