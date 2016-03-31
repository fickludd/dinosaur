package se.lth.immun

import java.io.File
import java.io.FileWriter
import java.io.BufferedWriter

object MzQuantML {

	val preAmble = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<MzQuantML xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
 xsi:schemaLocation="http://psidev.info/psi/pi/mzQuantML/1.0.0 ../../../schema/mzQuantML_1_0_0.xsd"
 xmlns="http://psidev.info/psi/pi/mzQuantML/1.0.0"
 version="1.0.0">
    <CvList>
        <Cv id="PSI-MS" uri="http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" version="2.25.0" fullName="Proteomics Standards Initiative Mass Spectrometry Vocabularies"/>
        <Cv id="UNIMOD" uri="http://www.unimod.org/obo/unimod.obo" fullName="Unimod"/>
    </CvList>
    <AnalysisSummary>
        <cvParam accession="MS:1001834" cvRef="PSI-MS" name="LC-MS label-free quantitation analysis"/>
        <cvParam accession="MS:1002019" cvRef="PSI-MS" value="true" name="label-free raw feature quantitation"/>
    </AnalysisSummary>
"""
	def inputFiles(inPath:String) = 
"""    <InputFiles>
        <RawFilesGroup id="the-file-group">
            <RawFile location="%s" id="the-file"/>
        </RawFilesGroup>
    </InputFiles>
""".format(inPath)

	def softwareList(version:String) = 
"""    <SoftwareList>
        <Software version="%s" id="Dinosaur"/>
    </SoftwareList>
    <DataProcessingList>
        <DataProcessing order="1" software_ref="Dinosaur" id="dino">
            <ProcessingMethod order="1"/>
        </DataProcessing>
    </DataProcessingList>
""".format(version)
	
	def midPart = 
"""    <AssayList id="AssayList_1">
        <Assay rawFilesGroup_ref="the-file-group" id="the-assay">
            <Label>
                <Modification>
                    <cvParam accession="MS:1002038" cvRef="PSI-MS" name="unlabeled sample"/>
                </Modification>
            </Label>
        </Assay>
    </AssayList>
	<FeatureList id="the-features" rawFilesGroup_ref="the-file-group">
"""
		
	def feature(ip:IsotopePattern, rtMap:Seq[Double], id:String) =
"""        <Feature charge="%d" mz="%.6f" rt="%.3f" id="%s">
            <MassTrace>%.3f %.6f %.3f %.6f</MassTrace>
        </Feature>
""".format(ip.z, ip.hills.head.total.centerMz, ip.apexHill.accurateApexRt(rtMap), id, ip.startRt(rtMap), ip.hills.head.total.centerMz, ip.endRt(rtMap), ip.hills.last.total.centerMz)
	
	val featureQuantPreAmble =
"""        <FeatureQuantLayer id="the-feature-intensities">
	        <ColumnDefinition>
	            <Column index="0">
	                <DataType>
	                    <cvParam cvRef="PSI-MS" accession="MS:1001141" name="intensity of precursor ion"/>
	                </DataType>
	            </Column>
		    </ColumnDefinition>
		    <DataMatrix>
"""
	def featureQuant(ip:IsotopePattern, id:String) = 
"""                <Row object_ref="%s">%f</Row>
""".format(id, ip.intensity)
		
	val postAmble = 
"""            </DataMatrix>
        </FeatureQuantLayer>
	</FeatureList>
<MzQuantML>"""
		
		
	
	def write(
			f:File, 
			features:Seq[IsotopePattern], 
			rtMap:Seq[Double],
			params:DinosaurParams
	):Long = {
		val t0 = System.currentTimeMillis
		
		val w = new BufferedWriter(new FileWriter(f))
		
		w.write(preAmble)
		w.write(inputFiles(params.mzML))
		w.write(softwareList(params.version))
		w.write(midPart)
		
		for (i <- 0 until features.length)
			w.write(feature(features(i), rtMap, "f"+i))
		
		w.write(featureQuantPreAmble)
			
		for (i <- 0 until features.length)
			w.write(featureQuant(features(i), "f"+i))
			
		w.write(postAmble)
		
		w.close
		System.currentTimeMillis() - t0
	}
}