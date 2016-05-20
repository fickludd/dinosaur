# Dinosaur

This repository contains the APACHE 2 licensed source code for the mass spectrometry MS1 feature detection tool Dinosaur. Dinosaur is an improved reimplementation of the feature finding algorithm of MaxQuant. For a description of the algorithm and how it was evaluated, we refer to the white paper

### Dinosaur: a refined open source peptide MS feature detector

###### Johan Teleman1,2, Aakash Chawade1, Marianne Sandin1, Fredrik Levander1,3, Johan Malmström2

1. Dept. of Immunotechnology, Lund University, Sweden
2. Dept. of Clinical Sciences Lund, Lund University, Sweden
3. Bioinformatics Infrastructure for Life Sciences (BILS), Lund University, Sweden

This README contains the following topics

- Getting started with Dinosaur
- The Dinosaur source code
- Dinosaur configuration and parameters
- Detailed algorithm explanation



Getting started with Dinosaur
=============================

To run Dinosaur `Java 1.6` or greater is needed. The only required input data is the MS data in the mzML-format. We recommend using ProteoWizard to produce the mzML. 

Dinosaur can be acquired either by downloading a pre-compiled release (click the Release tag), or by downloadning the source and compiling.

### Direct download of binary

Click the releases button and select the version you want (Click the .jar file).


### Get Dinosaur from source

Download the Dinosaur source using either the program `git`, or by direct download of a zip.

For git, install git and run the following on the commandline

> git clone https://github.com/fickludd/dinosaur

or navigate to https://github.com/fickludd/dinosaur and click 'Download ZIP' to download the source directly.



To execute Dinosaur
-------------------

Navigate to the Dinosaur jar, and run on the commandline

> java -jar Dinosaur-1.1.0.free.jar

where 1.1.0 should be replaced with the actual version you downloaded.

Running like this, without any parameters, will just output help for the available commands and how to run Dinosaur. To actually detect features on a mzML file, the following command would be suitable

> java -jar Dinosaur-1.1.0.free.jar --verbose --profiling --concurrency=4 myData.mzML

where myData.mzML should be replaced with your relevant file.



The Dinosaur source code
========================

Dinosaur is written in Scala 2.10.3 (http://scala-lang.org). This is a high level language that in compiled to produce programs for the java virtual machine (JVM). Compilation of Dinosaur and management of dependencies is achieved via maven. 

Program execution begins in the `Dinosaur` class, which handles the outermost program logic. Actual algorithm execution is thought performed in `FeatureFinder.analyzeMzML()`. Dinosaur main parameters are defined in `DinosaurParams`, and advanced algorithm parameters in `DinosaurAdvParams`. A DinosaurParams object is instantiated at the start of execution in the Dinosaur class, and passed around where ever parameters are needed.

Classes named SomethingReport hold logic related to plots in the plot trail. 

Compilation
-----------

Dinosaur requires several dependencies, as specified in the `pom.xml` file. These dependencies are archived in Maven central, and should be automagically downloaded by Maven.

| Dependency                                                          | Repository                       |
| ------------------------------------------------------------------- | -------------------------------- |
| MsNumpress                                                          | github.com/msnumpress/msnumpress | 
| Params, CLIApp                                                      | github.com/fickludd/treacle      | 
| Collections, Graphs, MsFeatureIO, MSFeatureProtocol, MzML, Proteins | github.com/fickludd/proteomicore | 


Dinosaur is compiled using `Maven 3.0.3` (https://maven.apache.org). Compile Dinosaur by running `mvn install` in the Dinosaur source code directory, which will generate a selfcontained executable and store it in the target folder, `target/Dinosaur-VERSION.free.jar`.


Dinosaur configuration and parameters
=====================================

All parameters, their description and default values can be listed directly from the executable jar. Running Dinosaur without parameters will present the most used parameters, which are related to output formats, meta-configuration, meta-data to report, plot trail configuration and such. By adding the flag --advHelp and some mzML file, advanced parameters relating to detailed algorithm execution will be listed and described.


  
Detailed algorithm explanation
==============================

Dinosaur performs feature detection in 4 steps. In the first step, profile spectra are centroided. In the second step, the centroid peaks of subsequent spectra are joined into hills, based on having sufficiently close m/z. In the third step, graphs are created by joining hills with a m/z difference close to ISOTOPE_MASS_DIFF / z and reasonable correlating intensity profiles. In the fourth step. The clusters from step 3 and filtered for charga-state consistency and isotope envelope intensities close to the corresponding averagine.

Step 1 - centroiding
--------------------

Centroiding is performed by detecting local maxima and extending these on both sides until a local minimum is observed (`Centroider.centroidSpectrum()`). For each such local peak, the m/z is computed as average of a gaussian fitted on the 3 most intense peaks (`Estimation.center()` and `Estimation.est3()`). The `maxIntensity` and `massEstPoints` advanced parameters can be used to tweak centroiding behavior.

Step 2 - hill building
----------------------

As spectra are centroided they are passed to the `HillBuilderActorParallel` for hill building. To speed up computation, this actor groups spectra into batches of 100 concurrent spectra (tunable by advanced parameter `hillBatchSize`), and build hills for each batch in parallell and soon as the batch is complete and there is available concurrency.

Batch hill building is performed by `RawHillActor`s. Here each new spectrum is added to the current set of hills via the `RawHillActor.processSpec()` and `RawHillActor.merge()` methods. The active hills are sorted by m/z and so is the centroid spectrum. The method iterates over both these lists simultaneously and compares the m/z difference of neighbouring hills and centroid peaks to a ppm threshold (advanced param `hillPPM`). If the m/z different is below the threshold the pair is added to a temporar buffer set. If the current m/z difference is outside the threshold, the current buffer is processed by iteratively removing the pair with the smallest m/z difference until no pairs are below the hillPPM threshold. For pairs that are remove, the centroid peak is appended to the respective hill profile. Centroid peaks that do not match to any hill are appended as new hills. 

After the merge operation, hills that have more that have their last centroid peak more than 1 (adv. `hillMaxMissing`) spectrum ago are moved from the active hill list to a completed hill list.

When hills are build for all spectrum batches, the resultant completed hill lists are sewn together so as to not split hills that happen to elute over spectrum batch borders. This process in highly similar to the `RawHillActor.merge()` method, but located in the `HillBuilderActorParallel.merge()`. 

For each hill in the sewn hill list, a smoothed intensity profile is calculated by a 3-point sliding median filter, followed by a 3-point sliding mean filter (adv. `hillSmoothMedianWindow` and `hillSmoothMeanWindow`). Following this, hills are recursive split in two in the rt dimension is a local minimum can be detected such that the smaller of the maximum value before the minimum and after the minimum is at least 30% more intense than the minimum (adv. `hillValleyFactor`). 

After hill splitting, hills of length 2 or 1 are discarded (adv. `hillMinLength`). For remaining hills, the total intensity is computed as the sum of the raw intensities, and the average m/z and m/z error are computed by bootstrap sampling 150 times (adv. `hillNBoots`) from the centroid peaks of the hill and computing weighted average m/z values. The average m/z and m/z error are computed from these bootstraps. 

Finally, to remove hills relating to ambient molecules, hill longer than 40 spectra (adv. `hillPeakFactorMinLength`) and with the first and last centroid peak larger than 50% of the hill apex are removed (adv. `hillPeakFactor`).

Step 3 - hill clustering
------------------------

Clustering hills is a somewhat complex operation if done in a speed efficient manner. The exact algorithm used is implemented in `Deisotoper.findEdges()` and `Deisotoper.ripEdges()`. In short, hills are sorted by end rt, and compared to one another in a limited fashion using a trailing index. A new datastructure in created with one list per hill, containing indices to other hills that it can be directly linked to. To collect the clusters, this structure is traversed while removing indices and adding them to a set. When all indices in the current set have been traversed the set cluster is complete. The cluster is finally represented by the set of detected edges between hills.

Two hills are joined by an edge if their m/z difference is close enough to the average isotope pattern mass difference divided by the considered charge (see `Deisotoper.fitsDiff()`), and if their intesity profiles have a cosine correlation >= 0.6 (adv. `deisoCorr`). Whether the m/z difference in close enought to the expected value is determined by the bootstrap-computed m/z errors of the hills and the error of having one sulphur isotope (S34) instead of two C13 isotopes.

Step 4 - cluster deconvolution
------------------------------

Like the spectra during hill building, clusters are grouped into groups of 1000 and deconvoluted in parallell. 

In each cluster, an iterative approach is employed where the 'best' found isotope pattern is removed until no more valid envelopes can be detected (`Cluster.deconvolve()`). The best current pattern is detected by seeding from the 100 most intense hills and each considered charge state. For each seed and charge state, hills are sought for at +- isotope pattern diff / z (`Cluster.extendSeed()`). If such hills are detected with a intesity profile cosine correlation >= 0.6 (adv. `deisoCorr`), the process is repeated adding another isotope pattern difference, until to matching hill is found. The detected primary isotope pattern is checked for local minima ( applying the `deisoValleyFactor` threshold ), which are truncated around the seed. For hills with masses < 1000, all extending hills that are of lower mass are removed. 

The final quality control of isotope patterns is the comparison vs. the isotope pattern of a averagine peptide of the same mass. The current isotope pattern is represented by a list of one intensity per hill. This list is compared to the intensities of the averagine isotopes down to 1e-6 prevalence by shifting the averagine isotopes and computing cosine correlations and percent explained averagine prevalence (`DinoUtil.optimalAlignOverlap()`). The explained prevalence and correlation are multiplied together and averagine shift with the highest value selected to determine the true monoisotopic peak. Finally the pattern is rejected if the averarine correlation is < 0.6 (adv. `averagineCorr`) or the explained prevalence < 0.5 (averagineExplained). 

  

