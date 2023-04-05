#! /usr/bin/env python

#  Copyright (C) 2015-16 Giovanni Scala, Istituto Nazionale di Fisica Nucleare (INFN),
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
import argparse
from methylUtils import preprocess


def probability(string):
    msg = "%r should be a real between 0 and 1" % string
    try:
        value = float(string)
        if 0.0 < value < 1.0:
            return value
        else:
            raise argparse.ArgumentTypeError(msg)
    except ValueError:
        raise argparse.ArgumentTypeError(msg)


def existentFile(string):
    if os.path.isfile(string):
        return string
    else:
        msg = "could not find file %s" % string
        raise argparse.ArgumentTypeError(msg)


parser = argparse.ArgumentParser()

parser.add_argument('--sample', action='store', dest='sample', required=True,
                    help='sample ID')
parser.add_argument('--primFile', action='store', dest='primersPath', required=True, type=existentFile,
                    help='file containing sequence IDs stored as "seqId,sequence,seqLen" one per row')
parser.add_argument('--readsFile', action='store', dest='fasta_path', required=True, type=existentFile,
                    help='fasta file containing reads to analyse')
parser.add_argument('--threshLen', action='store', dest='threshLen', type=probability, default=0.5,
                    help='the threshold constraint on the expected length (default=0.5)')
parser.add_argument('--primThresh', action='store', dest='primThresh', type=probability, default=0.8,
                    help='the percentage over the primer total number of bases of allowed mismatches for each primer (default=0.8)')

args = parser.parse_args()

sample = args.sample
primers_path = args.primersPath
fasta_path = args.fasta_path
threshLen = args.threshLen
primThresh = args.primThresh

preprocess(sample, primers_path, fasta_path, threshLen, primThresh)
