#  Copyright (C) 2015-16 Giovanni Scala, University of Naples Federico II,
#  Naples
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import multiprocessing
from methylUtils import extract_profiles, preprocess, make_summary, qualitative_Analysis
import argparse


def createQconf():
    """Yield successive n-sized chunks from l."""
    f_out = open("Qconf", "w")
    f_out.write("#! /bin/bash"+"\n"+"true")
    f_out.close()

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def probability(string):
    msg = "%r should be a real between 0 and 1" % string
    try:
        value = float(string)
        if 0.0 <= value <= 1.0:
            return value
        else:
            raise argparse.ArgumentTypeError(msg)
    except ValueError:
        raise argparse.ArgumentTypeError(msg)


def check_negative(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue


def existentFile(string):
    if os.path.isfile(string):
        return string
    else:
        msg = "could not find file %s" % string
        raise argparse.ArgumentTypeError(msg)


def existentDir(string):
    if os.path.isdir(string):
        if string.endswith("/"):
            return string
        else:
            return string + "/"
    else:
        msg = "could not find dir %s" % string
        raise argparse.ArgumentTypeError(msg)


def nonExistentDir(string):
    if not os.path.isdir(string):
        if string.endswith("/"):
            return string
        else:
            return string + "/"
    else:
        msg = "output folder already exists %s" % string
        raise argparse.ArgumentTypeError(msg)


# extracts methylation profiles from fasta file copntaing reads from sequencer in fasta format
def profiles_from_fasta(sample, primersPath, fasta_path, threshLen, primThresh, refFile, bisu_thresh, cgAmbThresh,
                        alignPropThresh, blastExecPath,
                        nThreadBlast, reward, penalty, soft_masking, dust, gapopen, gapextend, word_size, window_size):
    os.mkdir(sample)
    os.chdir(sample)
    preprocess(sample, primersPath, fasta_path, threshLen, primThresh)
    extract_profiles(refFile, primersPath, sample, bisu_thresh, bisuBaseThresh, cgAmbThresh, alignPropThresh,
                     blastExecPath,
                     nThreadBlast, reward, penalty, soft_masking, dust, gapopen, gapextend, word_size, window_size)
    os.chdir("..")

createQconf()

parser = argparse.ArgumentParser()

parser.add_argument('--refFile', action='store', dest='refFile', required=True, type=existentFile,
                    help='reference fasta file')

parser.add_argument('--primFile', action='store', dest='primersPath', required=True, type=existentFile,
                    help='file containing sequence IDs stored as "seqId,sequence,seqLen" one per row')

parser.add_argument('--bisu_thresh', action='store', dest='bisu_thresh', type=probability, default=0.98,
                    help='bisulfite efficiency threshold (min percentage of C (not CG) actually converted to T ) (default=0.98)')

parser.add_argument('--bisu_valid_Cs', action='store', dest='bisuBase', type=probability, default=0.0,
                    help='minimum percentage of Cs (not CpG) on which bisulfite efficiency is evaluated (default=0.0)')

parser.add_argument('--cgAmbThresh', action='store', dest='cgAmbThresh', type=probability, default=1.0,
                    help='maximum percentage of ambiguously aligned CG sites (C nor T are found on the sites (default=1.0)')

parser.add_argument('--alignPropThresh', action='store', dest='alignPropThresh', type=probability, default=0.6,
                    help='minimum proportion of  not clipped bases on the aligned read (default=0.6)')

parser.add_argument('--blastExecPath', action='store', dest='blastExecPath', type=existentDir, required=True,
                    help='path to the local blastn binary directory')

parser.add_argument('--nThread', action='store', dest='nThread', type=check_negative, default=1,
                    help='number of thread for execution (default=1)')

parser.add_argument('--reward', action='store', dest='reward', default=2,
                    help='blastn reward for a nucleotide match (see blastn documentation) (default=2)')

parser.add_argument('--penalty', action='store', dest='penalty', default=-3,
                    help='blastn penalty for a nucleotide mismatch (see blastn documentation) (default=-3)')
##soft_masking
parser.add_argument('--soft_masking', action='store', dest='soft_masking', default=True,
                    help='blastn soft_masking (default=True)')
#dust
parser.add_argument('--dust', action='store', dest='dust', default='yes',
                    help='blastn dust (default=20 64 1)')
parser.add_argument('--gapopen', action='store', dest='gapopen', default=5,
                    help='blastn cost to open a gap (default=5)')

parser.add_argument('--gapextend', action='store', dest='gapextend', default=2,
                    help='blastn cost to extend a gap (see blastn documentation) (default=2)')

parser.add_argument('--word_size', action='store', dest='word_size', default=11,
                    help='blastn ord size for initial match (see blastn documentation) (default=11)')

parser.add_argument('--window_size', action='store', dest='window_size', default=60,
                    help='blastn multiple hits window size (see blastn documentation) (default=60)')

parser.add_argument('--fastaDir', action='store', dest='fasta_path', required=True, type=existentDir,
                    help='folder containing fasta files with reads to analyse')

parser.add_argument('--threshLen', action='store', dest='threshLen', type=probability, default=0.5,
                    help='the threshold constraint on the expected length (default=0.5)')

parser.add_argument('--primThresh', action='store', dest='primThresh', type=probability, default=0.8,
                    help='the percentage over the primer total number of bases of allowed mismatches for each primer (default=0.8)')

parser.add_argument('--perfQual', action='store_true', dest='qualAn',
                    help='perform qualitative analisys. Requires qiime or macqiime installed on your system')
parser.add_argument('--qualDir', action='store', dest='qualDir', type=nonExistentDir,
                    help='qualitive analisys output folder')
parser.add_argument('--metaFile', action='store', dest='metaFile', type=existentFile,
                    help='sample metada file')
parser.add_argument('--qiime_path', action='store', dest='qiime_path', default="Qconf",type=existentFile,
                    help='qiime environmet script path. Not required if using anaconda installation')

args = parser.parse_args()

n_t = args.nThread
fastaDir = args.fasta_path
refFile = args.refFile
primers_path = args.primersPath

bisu_thresh = args.bisu_thresh
bisuBaseThresh = args.bisuBase
cgAmbThresh = args.cgAmbThresh
threshLen = args.threshLen
primThresh = args.primThresh

alignPropThresh = args.alignPropThresh
blastExecPath = args.blastExecPath
reward = args.reward
penalty = args.penalty
##soft_masking
soft_masking= args.soft_masking
dust=args.dust
gapopen = args.gapopen
gapextend = args.gapextend
word_size = args.word_size
window_size = args.window_size

qualAn = args.qualAn
summary = args.qualDir
qiime_env_path = args.qiime_path
sampleMeta = args.metaFile
fastaDir = args.fasta_path
# get list of samples names from fatsa files in fastaDir
samples = [each.split('.')[0] for each in os.listdir(fastaDir) if each.endswith('.fasta')]
#samples
# spawn computations over multiple CPUs

# determine maximum number of CPUs for each blast alignment
if n_t > len(samples):
    nThreadBlast = n_t.__floordiv__(len(samples))
else:
    nThreadBlast = 1

# process samples in chunks based on available resources
for sms in chunks(samples, n_t):
    jobs = []
    for sam in sms:
        fasta_path = fastaDir + "/" + sam + ".fasta"
        j = multiprocessing.Process(target=profiles_from_fasta, args=(
        sam, primers_path, fasta_path, threshLen, primThresh, refFile, bisu_thresh, cgAmbThresh, alignPropThresh,
        blastExecPath,
        nThreadBlast, reward, penalty, soft_masking, dust, gapopen, gapextend, word_size, window_size))
        jobs.append(j)
        j.start()

    for j in jobs:
        j.join()

# create global methylation profile summary by linking methylation profiles of all analysed samples

os.mkdir(summary)

# gather IDs of analysed regions
f = open(primers_path)
regions = set()
for line in f:
    fields = line.rstrip().split(",")
    regions.add(fields[0])
f.close()

make_summary(samples, regions, summary)

# if requested perform some qualitative analyses on computed methylation profiles
if qualAn:
    for reg in regions:
        destDir = summary + "/" + reg
        biom_reg = destDir + "/" + reg + ".biom"
        qualitative_Analysis(biom_reg, destDir, reg, sampleMeta, qiime_env_path)
