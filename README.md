# AmpliMethProfiler_RepetitiveElements
AmpliMethProfiler is a bioinformatic tool developed by Scala et al. (10.1186/s12859-016-1380-3) to perform average DNA methylation and epihaplotype analysis on deep targeted-bisulfite sequencing data.
For installation and usage guidelines of base version of AmpliMethProfiler tool please refer to repository indicated by Scala et al. at: https://sourceforge.net/projects/amplimethprofiler/.

Our group optimized the AmpliMethProfiler tool to perform DNA methylation analysis of Repetitive Elements (RE) including LINE-1, Alu and ribosomal DNA repeats using an Illumina-based, targeted-deep bisulfite sequencing pipeline. 
AmpliMethProfiler's DNA methylation calling is based on the alignment of sequencing reads on a specific reference sequence and the identification of mismatches at the position corresponding to the citosyne in a CpG dinucleotyde motif. Alignment is performed using blastn tool (please refer to Scala et al. manuscript for blastn standard settings) which integrates by default a low-complexity features masking module based on an heuristic algorythm called DUST, described by Morgulis et al. (10.1089/cmb.2006.13.1028)

Bisulfite conversion and subsequent target PCR-amplification of genomic DNA results in the conversion of all unmethylated cytosines (C) into thymines (T). 
For regions which are particularly rich in unmethylated Cs (both within or outside CpG contex). In genomic regions with high density of individual Cs and unmethylated CpG dinucleotides, bisulfite conversion results in a net reduction of sequence complexity.

We noticed that blastn default settings tagged and filtered out bisulfite-converted, unmethylated reads as low-complexity sequences. This resulted in a biased assessment of DNAm profile towards a highly methylated one.

We modified Scala's scripts by integrating an additional "--dust" argument to AmpliMethProfiler base command line which calls bioblast _dust_ property (https://biopython.org/docs/1.76/api/Bio.Blast.Applications.html?highlight=blastn#Bio.Blast.Applications.NcbiblastnCommandline). 
Dust argument allows a "yes" (default) or "no". Use the "no" argument to disable its usage.

**Usage:** 

1) Download the modified AmpliMethProfiler scripts from:
https://github.com/LabBrainAgeing/AmpliMethProfiler_RepetitiveElements/tree/main/AmpliMethProfiler.

2) Set up AmpliMethProfiler environment for your anaconda installation.
Guidelines to set up AmpliMethProfiler environment are provided at : https://sourceforge.net/projects/amplimethprofiler/

Alternatively, AmplyMethProfiler environment can be downloaded and installed from anaconda.org
  1) Using web interface, go to https://anaconda.org/ravaioli.cesco/AmpliMethProfilerEnv
  2) Select the environment, go to Files tab and click the file to download it
  3) Install the environment using the terminal:
  
	> conda env create $PWD/AmpliMethProfilerEnv.yml
	> conda activate AmpliMethProfilerEnv

3) Run base AmpliMethProfilerCommand line with DUST module disabled.
An example is provided at: https://github.com/LabBrainAgeing/AmpliMethProfiler_RepetitiveElements/blob/main/testData/Amplimeth_command.sh
