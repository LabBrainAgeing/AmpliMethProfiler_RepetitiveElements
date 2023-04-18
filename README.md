# AmpliMethProfiler_RepetitiveElements
AmpliMethProfiler is a bioinformatic tool developed by Scala et al. (10.1186/s12859-016-1380-3) to perform average DNA methylation and epihaplotype analysis on deep targeted-bisulfite sequencing data.
For installation and usage guidelines of base version of AmpliMethProfiler tool please follow the link to https://sourceforge.net/projects/amplimethprofiler/.

Our group optimized the AmpliMethProfiler tool to perform DNA methylation analysis of Repetitive Elements (RE) including LINE-1, Alu and ribosomal DNA repeats using an Illumina-based, targeted-deep bisulfite sequencing pipeline. 
AmpliMethProfiler's DNA methylation calling is based on the alignment of sequencing reads on a specific reference sequence and the identification of mismatches at the position corresponding to the citosyne in a CpG dinucleotyde motif. Alignment is performed using blastn tool (please refer to Scala et al. manuscript for blastn standard settings) which integrates by default a low-complexity features masking module based on an heuristic algorythm called DUST, described in Morgulis et al. (10.1089/cmb.2006.13.1028)

Bisulfite conversion and subsequent target PCR-amplification of genomic DNA results in the conversion of all unmethylated cytosines (C) into thymines (T). 
For regions which are particularly rich in unmethylated Cs (both within or outside CpG contex). In genomic regions with high density of individual Cs and unmethylated CpG dinucleotides, bisulfite conversion results in a net reduction of sequence complexity.

We noticed that blastn default settings filtered out these bisulfite-converted, generally unmethylated reads as low-complexity sequences resulting in a biased assessment towards higher levels of DNAm.

We modified Scala's scripts by integrating an additional "--DUST" argument to its base command line which calls bioblast _dust_ property (https://biopython.org/docs/1.76/api/Bio.Blast.Applications.html?highlight=blastn#Bio.Blast.Applications.NcbiblastnCommandline). DUST argument allows a "yes" (default) or "no" format.


**Usage:** 

