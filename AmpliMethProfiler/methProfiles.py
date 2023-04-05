#! /usr/bin/env python

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

from methylUtils import extract_profiles
import os
import argparse


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


def existentDir(string):
    if os.path.isdir(string):
        if string.endswith("/"):
            return string
        else:
            return string + "/"
    else:
        msg = "could not find dir %s" % string
        raise argparse.ArgumentTypeError(msg)


def existentFile(string):
    if os.path.isfile(string):
        return string
    else:
        msg = "could not find file %s" % string
        raise argparse.ArgumentTypeError(msg)


parser = argparse.ArgumentParser()

parser.add_argument('--refFile', action='store', dest='refFile', required=True, type=existentFile,
                    help='reference fasta file')
parser.add_argument('--primFile', action='store', dest='primersPath', required=True, type=existentFile,
                    help='file containing sequence IDs stored as "seqId,sequence,seqLen" one per row')
parser.add_argument('--sample', action='store', dest='sample', required=True,
                    help='sample name')
parser.add_argument('--bisu_thresh', action='store', dest='bisu_thresh', type=probability, default=0.98,
                    help='bisulfite efficiency threshold (min percentage of C (not CpG) actually converted to T) (default=0.98)')
parser.add_argument('--bisu_valid_Cs', action='store', dest='bisuBase', type=probability, default=0.0,
                    help='minimum percentage of Cs (not CpG) on which bisulfite efficiency is evaluated (default=0.0)')
parser.add_argument('--cgAmbThresh', action='store', dest='cgAmbThresh', type=probability, default=1.0,
                    help='maximum percentage of ambiguously aligned CG sites (C nor T are found on the sites (default=1.0)')
parser.add_argument('--alignPropThresh', action='store', dest='alignPropThresh', type=probability, default=0.6,
                    help='minimum proportion of  not clipped bases on the aligned read (default=0.6)')
parser.add_argument('--blastExecPath', action='store', dest='blastExecPath', required=True, type=existentDir,
                    help='path to the local blastn binary directory')
parser.add_argument('--nThread', action='store', dest='nThreadBlast', type=check_negative, default=1,
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

args = parser.parse_args()

refFile = args.refFile
primersPath = args.primersPath
sample = args.sample
bisu_thresh = args.bisu_thresh
bisuBaseThresh = args.bisuBase
cgAmbThresh = args.cgAmbThresh
alignPropThresh = args.alignPropThresh
blastExecPath = args.blastExecPath
nThreadBlast = args.nThreadBlast
reward = args.reward
penalty = args.penalty
soft_masking = args.soft_masking
dust=args.dust
gapopen = args.gapopen
gapextend = args.gapextend
word_size = args.word_size
window_size = args.window_size

extract_profiles(refFile, primersPath, sample, bisu_thresh, bisuBaseThresh, cgAmbThresh, alignPropThresh, blastExecPath,
                 nThreadBlast,
                 reward, penalty, soft_masking, dust, gapopen, gapextend, word_size, window_size)
