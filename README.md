# AmpliMethProfiler_RepetitiveElements
AmpliMethProfiler is a bioinformatic tool developed by Scala et al. (10.1186/s12859-016-1380-3) to perform average DNA methylation and epihaplotype analysis on deep targeted-bisulfite sequencing data.
For installation and usage of base version of AmpliMethProfiler tool please refer to https://sourceforge.net/projects/amplimethprofiler/.

Our group optimized the AmpliMethProfiler tool to perform DNA methylation analysis of DNA Repetitive Elements including LINE-1, Alu and ribosomal DNA repeats. 

We modified Scala's scripts by integrating an additional "--DUST" argument to its base command line. DUST is a BLAST module used to mask low-complexity sequences of nucleotide queries (10.1089/cmb.2006.13.1028). 

This allows to limit read loss as bisulfite-conversion reduces sequence complexity by converting unmethylated Cytosines (C) into Thymines (T).

Usage: 
