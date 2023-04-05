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

import argparse, os
from methylUtils import qualitative_Analysis

def createQconf():
    """Yield successive n-sized chunks from l."""
    f_out = open("Qconf", "w")
    f_out.write("#! /bin/bash"+"\n"+"true")
    f_out.close()

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


parser = argparse.ArgumentParser()
parser.add_argument('--outDir', action='store', dest='qualDir', type=nonExistentDir,
                    help='analisys output folder')
parser.add_argument('--metaFile', action='store', dest='metaFile', type=existentFile,
                    help='sample metadata file')
parser.add_argument('--qiime_path', action='store', dest='qiime_path', default="Qconf", type=existentFile,
                    help='qiime environmet script path. Not required if using anaconda installation')
parser.add_argument('--biom_File', action='store', dest='biom_reg', type=existentFile,
                    help='biom file containing methylation profile counts for each sample')
parser.add_argument('--regName', action='store', dest='region_name',
                    help='id of the analysed region')

args = parser.parse_args()

qiime_env_path = args.qiime_path
sampleMeta = args.metaFile
region = args.region_name
destDir = args.qualDir
biom_reg = args.biom_reg
os.mkdir(destDir)
# performs a series of qualitative analises using local qiime installation for each one of the analyzed regions
qualitative_Analysis(biom_reg, destDir, region, sampleMeta, qiime_env_path)
