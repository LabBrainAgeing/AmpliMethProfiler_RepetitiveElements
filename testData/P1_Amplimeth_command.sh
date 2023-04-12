#!/bin/sh
mkdir ~/Desktop/RiboDown_NoDust/P1/P1_AMPLIMETH/FlunkrDNA1
cd ~/Desktop/RiboDown_NoDust/P1/P1_AMPLIMETH/FlunkrDNA1/
python ~/Dropbox/Sequencing/Meth_Seq/EpiAlleles/AmpliMethProfiler/ampliMethProfiler.py --refFile ~/Desktop/RiboDown_NoDust/P1/P1_ReferenceSeq/FlunkrDNA1_ReferenceSeq.fasta --primFile ~/Desktop/RiboDown_NoDust/P1/P1_PrimFile/FlunkrDNA1_PrimFile.csv --blastExecPath /Users/MGiulia/anaconda2/envs/AmpliMethProfilerEnv/bin --fastaDir ~/Desktop/RiboDown_NoDust/P1/P1_FASTA --perfQual --dust no --metaFile ~/Desktop/RiboDown_NoDust/P1/Samplesheet_RiboDown_P1.txt --qualDir ~/Desktop/RiboDown_NoDust/P1/P1_AMPLIMETH/FlunkrDNA1/Stat --threshLen 0.2 --nThread 10
mkdir ~/Desktop/RiboDown_NoDust/P1/P1_AMPLIMETH/FlunkrDNA2
cd ~/Desktop/RiboDown_NoDust/P1/P1_AMPLIMETH/FlunkrDNA2/
python ~/Dropbox/Sequencing/Meth_Seq/EpiAlleles/AmpliMethProfiler/ampliMethProfiler.py --refFile ~/Desktop/RiboDown_NoDust/P1/P1_ReferenceSeq/FlunkrDNA2_ReferenceSeq.fasta --primFile ~/Desktop/RiboDown_NoDust/P1/P1_PrimFile/FlunkrDNA2_PrimFile.csv --blastExecPath /Users/MGiulia/anaconda2/envs/AmpliMethProfilerEnv/bin --fastaDir ~/Desktop/RiboDown_NoDust/P1/P1_FASTA --perfQual --dust no --metaFile ~/Desktop/RiboDown_NoDust/P1/Samplesheet_RiboDown_P1.txt --qualDir ~/Desktop/RiboDown_NoDust/P1/P1_AMPLIMETH/FlunkrDNA2/Stat --threshLen 0.2 --nThread 10
