#!/bin/sh
mkdir $PWD/RiboProm1
cd $PWD/RiboProm1
python $PWD/ampliMethProfiler.py --refFile $PWD/RiboProm1_ReferenceSeq.fasta --primFile $PWD/RiboProm1_PrimFile.csv --blastExecPath /home/user/anaconda2/envs/AmpliMethProfilerEnv/bin --fastaDir $PWD/fasta --perfQual --dust no --metaFile $PWD/Samplesheet.txt --qualDir $PWD/RiboProm1/Stat --threshLen 0.2 --nThread 10
