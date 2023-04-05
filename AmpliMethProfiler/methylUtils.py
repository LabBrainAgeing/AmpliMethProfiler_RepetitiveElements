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
import re
import ntpath
import difflib
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
import time


# recognise the reference sequence by comparing the sequence in molDict with the 5' and 3' ends of the read
# search for matches with the given primers in the 5' side of the sequence allowing for mismatches
# the percentage of allowed mismatches for each primer is given as a threshold value between 1 and primThresh
# It starts searching with threshold 1 (perfect sequence match) and decrease the threshold value until a match is found
# or the maximum number of allowed mismatches is reached
def find_targSeq(line, primThresh, molDict):
    primList = []

    thresh5 = 1
    found5 = False

    primLens = list(set([len(a) for a in molDict.keys()]))
    nPrims = len(primLens)

    while thresh5 >= primThresh and not found5:
        i = 0
        pk = []
        # for each given primer try to match the first part of the sequence
        while i < nPrims:
            primerFound = line[0:primLens[i]]
            primerKey = difflib.get_close_matches(primerFound, molDict.keys(), 1, thresh5)
            if primerKey:
                pk.append(primerKey)
            i = i + 1

        # check if one or more matching primers were found
        if (pk):
            found5 = True
            primList5 = []
            # record found primers
            for k in pk:
                primerKey = k[0]
                primList5.append(molDict[primerKey])
                primList5 = list(set(primList5))
        else:
            # if no primer matched the beginning of the sequence lower the threshold
            thresh5 = thresh5 - 0.1

    # same as 5' side
    thresh3 = 1
    found3 = False

    while thresh3 >= primThresh and not found3:
        i = 0
        pk = []
        while i < nPrims:
            primerFound = line[(len(line) - primLens[i] - 1):(len(line) - 1)]
            primerKey = difflib.get_close_matches(primerFound, molDict.keys(), 1, thresh3)
            if (primerKey):
                pk.append(primerKey)
            i = i + 1

        if pk:
            found3 = True
            primList3 = []
            for k in pk:
                primerKey = k[0]
                primList3.append(molDict[primerKey])
                primList3 = list(set(primList3))
        else:
            thresh3 = thresh3 - 0.1

    if found5 and found3:
        # if target sequences are found in both sides of the read then chose the ones with less mismatches
        primList = primList5 if (thresh5 >= thresh3) else primList3
    elif (found5):
        primList = primList5
    elif (found3):
        primList = primList3

    return primList


# check if the observed length of th input read meets the threshold constraint on the expected length
def filter_len(read_id, obsLen, threshLenRead, targSeqDict):
    # extract the expected targSeq length from the dictionary
    expLentargSeq = float(targSeqDict[read_id])
    # calculate the minimum length accepted
    minLen = expLentargSeq - (expLentargSeq * threshLenRead)
    # calculate the maximum length accepted
    maxLen = expLentargSeq + (expLentargSeq * threshLenRead)

    if (obsLen >= minLen) and (obsLen <= maxLen):
        seq = True
    else:
        seq = False
    return seq


# finds c (not cg) and cg position in an input gapped sequence
def find_c_cg_pos(ref_seq):
    ref_seq = ref_seq.upper()
    cg_pos = [m.start() for m in re.finditer('C-*G', ref_seq)]
    c_pos = [m.start() for m in re.finditer('C(?=-*[^-G])', ref_seq)]
    return cg_pos, c_pos


# remaps position from an ungapped sequence into the same gapped sequence
def remap_pos_align(ref_seq, old_pos):
    ref_seq = ref_seq.upper()
    gap_pos = [m.start() for m in re.finditer('-', ref_seq)]
    new_cpos = []
    gap_ind = 0
    for cp in old_pos:
        while gap_ind < len(gap_pos) and gap_pos[gap_ind] <= cp:
            gap_ind += 1
        new_cpos.append(cp + gap_ind)
    return new_cpos


# bisulfite covert the reference sequence corresponding to the given id (referenceSeqId)
# reads the reference sequences from an input fasta file and store the result in fasta format using the given output path


def bisuReferenceFasta(refPath, referenceSeqId, bisuPath):
    f = open(refPath)
    f_out = open(bisuPath, "w")

    head = f.readline()
    # search the targSeq sequence by looking for correspondence of the given id with the id in the fasta header
    while (len(head) > 0) & (head.rstrip().replace(">", "") != referenceSeqId):
        f.readline()
        head = f.readline()

    if head:
        f_out.write(head)
        seq = f.readline().rstrip().upper()
        # converts all C (not CG) bases in T
        bisu_str = re.sub("C(?!G)", "T", seq)
        f_out.write(bisu_str + "\n")

    f.close()
    f_out.close()


# align the given sequence file against a reference by using the local blast installation and xml as output format
# seq_filePath  ->  fasta file containing input query sequences
# ref_file      ->  fasta file containing the subject sequence (the reference)
# align_name    ->  alignment output name
# blastExecPath ->  local blast executables folder path
# n_thread      ->  number of thread to use for the BLAST execution
# ....          ->  other blastn parameters
def blastAlign(seq_filePath, ref_file, align_name, blastExecPath, n_thread,
               reward=2, penalty=-3, gapopen=5, gapextend=2, word_size=11, window_size=60, outfmt=5, soft_masking=True, dust="yes"):
    os.mkdir("blastDB")
    os.chdir("blastDB")

    # build the blast DB
    os.system(
        blastExecPath + "makeblastdb -in " + ref_file + " -title=" + align_name +
        " -dbtype nucl -parse_seqids -out ./" + align_name + " -logfile blastDB.log")
    os.chdir("../")

    # build the blastn command line
    blastn_cline = NcbiblastnCommandline(cmd=blastExecPath + 'blastn', query=seq_filePath,
                                         db="blastDB/" + "/" + align_name,
                                         max_target_seqs=1, reward=2, penalty=-3, gapopen=5, gapextend=2,
                                         word_size=11, window_size=60, outfmt=5, num_threads=n_thread, soft_masking=True, dust=dust,
                                         out=align_name + ".xml")
    stdout, stderr = blastn_cline()


# analyzes blast output and computes for passing filter reads:
#       - methylation profiles
#       - alignment in clear format
#       - summary/quality statistics
# INPUT:
# blast_out         -> blastn output in xml format
# referenceSeqId    -> reference sequence ID
# refFile           -> fasta file containing reference sequences
# out_name          -> analysis output name prefix
# bisu_eff_thresh   -> bisulfite efficiency threshold (min percentage of C (not CG) actually converted to T )
# bisuBaseThresh    -> minimum percentage of reference non-CpG cytosine residues to be assayed in order to consider the bisulfite efficiency estimate valid
# cg_amb_thresh     -> maximum percentage of ambiguously aligned CG sites (C nor T are found on the sites)
# minLenThresh      -> minimum proportion of  not clipped bases on the aligned read
def getProfilesFromBlast(blast_out, referenceSeqId, refFile, out_name, bisu_eff_thresh=0, cg_amb_thresh=1,
                         minLenThresh=0.6, bisuBaseThresh=0.0):
    # load a dictionary with c and cg pos for each region
    # load a dictionary with bisulphinated reference sequences
    targSeqSeqCGpos = {}
    targSeqSeq = {}

    f = open(refFile)
    line = f.readline().rstrip()
    while line:
        # retrieve region ID from the header
        seqNam = line.replace(">", "")
        # read the reference sequence
        seq = f.readline().rstrip().upper()
        # extract reference c and cg positions from the reference sequence
        targSeqSeqCGpos[seqNam] = find_c_cg_pos(seq)
        targSeqSeq[seqNam] = seq
        line = f.readline().rstrip()
    f.close()

    c_pos = targSeqSeqCGpos[referenceSeqId][1]
    cg_pos = targSeqSeqCGpos[referenceSeqId][0]

    # store percent of non deaminated bases for each c or cg site
    c_meth_tot = [0.0 for i in c_pos]
    cg_meth_tot = [0.0 for i in cg_pos]

    # store the number of assessed bases for each c or cg locus
    n_valid_sites_c = [0 for i in c_pos]
    n_valid_sites_cg = [0 for i in cg_pos]
    # store the number of passing filter reads
    n_valid_reads = 0
    # store the total number of analysed reads
    tot_reads = 0

    f_out = open(out_name, "w")
    f_out_align = open(out_name + ".align", "w")
    f_out_stats = open(out_name + ".stats", "w")

    # parse BLASTn output
    result = open(blast_out, "r")
    blast_records = NCBIXML.parse(result)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            hsp = alignment.hsps[0]

            tot_reads += 1
            aln_ref = hsp.sbjct
            aln_seq = hsp.query
            match = hsp.match

            # if necessary, reverse complemented the alignment
            if hsp.sbjct_start > hsp.sbjct_end:
                aln_ref = str(Seq(aln_ref).reverse_complement())
                aln_seq = str(Seq(aln_seq).reverse_complement())
                match = match[::-1]
                tmp = hsp.sbjct_start
                hsp.sbjct_start = hsp.sbjct_end
                hsp.sbjct_end = tmp

            # fill with N clipped bases from the query and the subject sequence
            start_gaps = "N" * (hsp.sbjct_start - 1)
            end_gaps = "N" * (len(targSeqSeq[referenceSeqId]) - hsp.sbjct_end)

            aln_ref = start_gaps + aln_ref + end_gaps
            aln_seq = start_gaps + aln_seq + end_gaps
            match = start_gaps + match + end_gaps

            # find c and cg position in the gapped sequences from the alignment
            ref_c_pos = remap_pos_align(aln_ref, c_pos)
            ref_cg_pos = remap_pos_align(aln_ref, cg_pos)

            # evaluate c and cg conversion profiles for the aligner read
            # 1 means non deaminated c ; 0 means C->T transition; 2 for ambiguity C->[A,G,-]
            cg_met = [1 if aln_seq[i] == 'C' else 0 if aln_seq[i] == 'T' else 2 for i in ref_cg_pos]
            c_met = [1 if aln_seq[i] == 'C' else 0 if aln_seq[i] == 'T' else 2 for i in ref_c_pos]

            # compute bisulfite efficiency for the aligned read as 1 minus the ratio of deaminated C
            if (c_met.count(1) + c_met.count(0)) > 0:
                bisu_eff = 1 - (float(c_met.count(1)) / float(c_met.count(1) + c_met.count(0)))
            else:
                bisu_eff = 0

            bisuBase = float(c_met.count(1) + c_met.count(0)) / float(len(c_met))

            # compute the ratio of ambiguous CG: C->[A,G,-] on the aligned read
            uncertainCG = float(cg_met.count(2)) / len(cg_met)

            # compute the proportion of not clipped bases on the aligned read
            alignProp = float(hsp.align_length) / len(targSeqSeq[referenceSeqId])

            # discard reads that do not satisfy input filtering criteria
            if bisu_eff >= bisu_eff_thresh and bisuBase >= bisuBaseThresh and uncertainCG <= cg_amb_thresh and alignProp >= minLenThresh:
                n_valid_reads += 1

                # sum up the number of non ambiguous sites
                n_valid_sites_c = [x + y for x, y in zip(n_valid_sites_c, [1 if c != 2 else 0 for c in c_met])]
                n_valid_sites_cg = [x + y for x, y in zip(n_valid_sites_cg, [1 if c != 2 else 0 for c in cg_met])]

                # set to zero ambiguous sites in the read
                c_met_sum = [c if c != 2 else 0.0 for c in c_met]
                cg_met_sum = [cg if cg != 2 else 0.0 for cg in cg_met]

                # sum up the number of non deaminated sites
                c_meth_tot = [x + y for x, y in zip(c_meth_tot, c_met_sum)]
                cg_meth_tot = [x + y for x, y in zip(cg_meth_tot, cg_met_sum)]

                # write  the methylation profile for the read
                f_out.write(" ".join(map(str, cg_met)) + "\n")

                # write  the alignment in clear format
                f_out_align.write(">" + blast_record.query + "\n")
                f_out_align.write(
                    "bisulfite efficiency: " + str(bisu_eff) + " % of assayed Cs: " + str(bisuBase) + "\n")
                f_out_align.write(aln_ref + "\n")
                f_out_align.write(match + "\n")
                f_out_align.write(aln_seq + "\n")
    # compute the total percent of not deaminated Cs an CGs over all the aligned reads
    c_meth_tot = [1 - (x / y) for x, y in zip(c_meth_tot, [nc if nc > 0 else -1 for nc in n_valid_sites_c])]
    cg_meth_tot = [(x / y) for x, y in zip(cg_meth_tot, [cg if cg > 0 else -1 for cg in n_valid_sites_cg])]

    # write some summary statistics on the analysed alignment
    f_out_stats.write(str(tot_reads) + " : Analyzed reads \n")
    f_out_stats.write(str(n_valid_reads) + " : Passing filter reads \n\n")

    f_out_stats.write("Meth percent \n")
    f_out_stats.write("position " + " ".join(map(str, cg_pos)) + "\n")
    f_out_stats.write(" ".join(map(str, cg_meth_tot)) + "\n")
    f_out_stats.write("# of valid Cs " + " ".join(map(str, n_valid_sites_cg)) + "\n\n")

    f_out_stats.write("Bisu efficiency per site \n")
    f_out_stats.write("position " + " ".join(map(str, c_pos)) + "\n")
    f_out_stats.write(" ".join(map(str, c_meth_tot)) + "\n")
    f_out_stats.write("# of valid Cs " + " ".join(map(str, n_valid_sites_c)) + "\n\n")

    f_out.close()
    f_out_stats.close()
    f_out_align.close()
    result.close()


# open output file handlers, one for each target region and store them in a dictionary using sample-regionId as key
def openSampletargSeqHandlers(seqList, sam):
    outFiles = {}
    for seq in seqList:
        outFiles[(seq, sam)] = (open(seq + "_" + sam + ".fasta", "w"), 0)
    return outFiles


# close all output file handlers
def closeSampletargSeqHandlers(file_dict):
    for file in file_dict.values():
        file[0].close()
    for file in file_dict.values():
        if file[1] == 0:
            os.delete(file.name)


# load a dictionary where primers sequences are associated with target sequence IDs
# read records from a file where info is stored as "seqId,sequence,seqLen" one per row
def preprocess(sample, primers_path, fasta_path, threshLen=0.5, primThresh=0.8):
    f = open(primers_path)
    molDict = {}
    for line in f:
        fields = line.rstrip().split(",")
        molDict[fields[1]] = fields[0]
    f.close()
    # load a dictionary where the targSeq name is associated with the target sequence length
    f = open(primers_path)
    targSeqLenDict = {}
    for line in f:
        fields = line.split(",")
        targSeqLenDict[fields[0]] = fields[2]
    f.close()
    # output file
    f_outs = openSampletargSeqHandlers(targSeqLenDict.keys(), sample)
    # read single line fasta file
    f = open(fasta_path)
    notargSeqPrimerMatched = 0
    multiplePrimerMatched = 0
    # read first header line
    headLine = f.readline().rstrip()
    while headLine:
        primerKey = ""
        header = headLine.split(" ")

        line = f.readline().rstrip()
        lineSeq = True

        while lineSeq:
            nextLine = f.readline().rstrip()
            if (not nextLine) or nextLine[0] == ">":
                lineSeq = False
                headLine = nextLine
            else:
                line += nextLine

        # recompute the length of the sequence
        newLen = len(line)

        primList = find_targSeq(line=line, primThresh=0.8, molDict=molDict)

        # no target sequence found with the given threshold
        if len(primList) == 0:
            notargSeqPrimerMatched += 1
        elif len(primList) > 1:
            multiplePrimerMatched += 1
        else:
            targSeq = primList[0]
            # if the read passes length filters then add it to the output
            if filter_len(targSeq, newLen, threshLen, targSeqDict=targSeqLenDict):
                # add information to the header
                header.append("length=" + str(newLen))
                header.append("sample=" + sample)
                header.append("targSeqPrimer=" + targSeq)
                f_outs[(targSeq, sample)][0].write(" ".join(header) + "\n")
                f_outs[(targSeq, sample)][0].write(line + "\n")
                f_outs[(targSeq, sample)] = (f_outs[(targSeq, sample)][0], f_outs[(targSeq, sample)][1] + 1)
    f.close()
    closeSampletargSeqHandlers(f_outs)


# extracts methylation profiles from fasta file using BLAST alignment
def extract_profiles(refFile, primersPath, sample, bisu_thresh, bisuBaseThresh, cgAmbThresh, alignPropThresh,
                     blastExecPath,
                     nThreadBlast, reward, penalty, soft_masking, dust, gapopen, gapextend, word_size, window_size):
    wdir = os.getcwd()
    # load sequence IDs from a file where info is stored as "seqId,PrimerSequence,seqLen" one per row
    f = open(primersPath)
    molDict = []
    for line in f:
        fields = line.rstrip().split(",")
        molDict.append(fields[0])
    f.close()
    # perform the analysis for each sequenced target region
    for targSeq in set(molDict):
        out_name = targSeq + "_" + sample

        seqFile = wdir + "/" + out_name + ".fasta"
        # store results in a new directory for each region-sample
        os.mkdir(out_name)
        os.chdir(out_name)
        # compute and save the bisulfite converted reference for the region
        bisuRef = wdir + "/" + out_name + "/" + targSeq + ".bisu"
        bisuReferenceFasta(refFile, targSeq, bisuRef)
        # align the input sequences against the bisulphinated reference
        blastAlign(seqFile, bisuRef, out_name, blastExecPath, nThreadBlast, reward=reward, penalty=penalty,
                   soft_masking=soft_masking, dust=dust, gapopen=gapopen,
                   gapextend=gapextend, word_size=word_size, window_size=window_size)

        blastOut = wdir + "/" + out_name + "/" + out_name + ".xml"

        profFile = out_name + ".out"
        # compute methylation profile and statistics for the region-sample
        getProfilesFromBlast(blastOut, targSeq, refFile, profFile, bisu_thresh, cgAmbThresh, alignPropThresh,
                             bisuBaseThresh)

        biomFile = out_name + ".biom"
        makeSingleBiom(profFile, biomFile, sample)

        profTab = out_name + ".csv"
        makeSingleAbTab(profFile, profTab, sample)

        os.chdir("../")


# converts output profiles in the sparse biom format ("http://biom-format.org") for metagenomic based analyses
# INPUT:
# profFile    -> .out file containing methylation profiles of the sample generated with getProfilesFromBlast.py module
# sample      -> sample name
# out_name    -> output file name
def makeSingleBiom(profFile, out_name, sample):
    profiles = {}
    f = open(profFile)
    biom_out = open(out_name, "w")

    # counts occurences of each methylation profile in the sample
    for line in f:
        line = line.rstrip()
        if not line.__contains__("2"):
            if line in profiles:
                profiles[line] += 1
            else:
                profiles[line] = 1

    keys = sorted(profiles.keys())
    nSpec = len(keys)

    # save methylation profiles in biom format, we exploit the taxonomy metadata field to classify profiles based on their number
    # of methylated CpG
    biom_out.write(
        '{"id":null,\n"format": "Biological Observation Matrix 0.9.1",\n"format_url": "http://biom-format.org",'
        '\n"type": "OTU table",\n"generated_by": "AmpliMethProfiler",\n"date": "' + str(time.strftime("%d/%m/%Y %X")) +
        '",\n"rows":[')

    for i in xrange(nSpec - 1):
        n_met = sum([int(val) for val in keys[i].replace(" ", "")])
        biom_out.write('\n{"id":"' + keys[i].replace(" ", "") + '", "metadata":{"taxonomy":["' + str(n_met) + '"]}},')

    n_met = sum([int(val) for val in keys[i].replace(" ", "")])
    biom_out.write(
        '\n{"id":"' + keys[nSpec - 1].replace(" ", "") + '", "metadata":{"taxonomy":["' + str(n_met) + '"]}}')
    biom_out.write("\n],")

    biom_out.write('\n"columns": [\n{"id":"' + sample + '", "metadata":null}\n],\n')
    biom_out.write('\n"matrix_type": "sparse",\n"matrix_element_type": "int",')
    biom_out.write('\n"shape": [' + str(nSpec) + ' ,1 ],\n"data":[\n')

    for i in xrange(nSpec - 1):
        biom_out.write('[' + str(i) + ',0,' + str(profiles[keys[i]]) + '],\n')
    biom_out.write('[' + str(nSpec - 1) + ',0,' + str(profiles[keys[nSpec - 1]]) + ']\n')
    biom_out.write("]\n}")
    biom_out.close()


# save methylation profiles in a profile abundances file, for each profile we record the number of methylated CpG and
# the occurences of that profile in the sample
# INPUT:
# profFile    -> .out file containing methylation profiles of the sample generated with getProfilesFromBlast.py module
# sample      -> sample name
# out_name    -> output file name
def makeSingleAbTab(profFile, out_name, sample):
    profiles = {}
    f = open(profFile)
    tab_out = open(out_name, "w")

    for line in f:
        line = line.rstrip()
        if not line.__contains__("2"):
            if line in profiles:
                profiles[line] += 1
            else:
                profiles[line] = 1

    keys = sorted(profiles.keys())
    nSpec = len(keys)

    tab_out.write('id\tprofile\tcount\tn_meths\n')

    for i in xrange(nSpec):
        n_met = sum([int(val) for val in keys[i].replace(" ", "")])
        tab_out.write(
            sample + '\t' + keys[i].replace(" ", "") + '\t' + str(profiles[keys[i]]) + '\t' + str(n_met) + '\n')

    tab_out.close()


# converts output profiles in the sparse biom format ("http://biom-format.org") for metagenomic based analyses
# INPUT:
# sfiles    -> list of .out file paths containing methylation profiles of each sample generated with getProfilesFromBlast.py module
# out_name  -> output file name
def makeBiom(sfiles, out_name):
    sample = []
    sample += [re.sub(".out", "", ntpath.basename(name).split('_', 1)[-1]) for name in sfiles]

    nSam = len(sample)

    biom_out = open(out_name, "w")
    profiles = {}
    for i in range(nSam):
        file = open(sfiles[i])
        # counts occurences of each methylation profile in the sample
        for line in file:
            line = line.rstrip()
            if not line.__contains__("2"):
                if line in profiles:
                    profiles[line][i] += 1
                else:
                    profiles[line] = [0] * nSam
                    profiles[line][i] = 1
    keys = sorted(profiles.keys())
    nSpec = len(keys)

    # save methylation profiles in biom format, we exploit the taxonomy metadata field to classify profiles based on their number
    # of methylated CpG
    biom_out.write(
        '{"id":null,\n"format": "Biological Observation Matrix 0.9.1",\n"format_url": "http://biom-format.org",'
        '\n"type": "OTU table",\n"generated_by": "AmpliMethProfiler",\n"date": "' + str(time.strftime("%d/%m/%Y %X")) +
        '",\n"rows":[')

    if nSpec > 1:
        for i in xrange(nSpec - 1):
            n_met = sum([int(val) for val in keys[i].replace(" ", "")])
            biom_out.write(
                '\n{"id":"' + keys[i].replace(" ", "") + '", "metadata":{"taxonomy":["' + str(n_met) + '"]}},')
    n_met = sum([int(val) for val in keys[nSpec - 1].replace(" ", "")])
    biom_out.write(
        '\n{"id":"' + keys[nSpec - 1].replace(" ", "") + '", "metadata":{"taxonomy":["' + str(n_met) + '"]}}')
    biom_out.write("\n],")

    biom_out.write('\n"columns": [\n')

    if nSam > 1:
        for i in xrange(nSam - 1):
            biom_out.write('{"id":"' + sample[i] + '", "metadata":null},\n')
    biom_out.write('{"id":"' + sample[len(sample) - 1] + '", "metadata":null}')

    biom_out.write('\n],\n"matrix_type": "sparse",\n"matrix_element_type": "int",')
    biom_out.write('\n"shape": [' + str(nSpec) + ' ,' + str(nSam) + '],\n"data":[\n')

    for i in xrange(nSpec - 1):
        for j in xrange(nSam):
            biom_out.write('[' + str(i) + ',' + str(j) + ',' + str(profiles[keys[i]][j]) + '],\n')
    if nSam > 1:
        for j in xrange(nSam - 1):
            biom_out.write('[' + str(nSpec - 1) + ',' + str(j) + ',' + str(profiles[keys[nSpec - 1]][j]) + '],\n')
    biom_out.write('[' + str(nSpec - 1) + ',' + str(nSam - 1) + ',' + str(profiles[keys[nSpec - 1]][nSam - 1]) + ']\n')

    biom_out.write("]\n}")
    biom_out.close()


# save methylation profiles in a profile abundances file, for each profile we record the number of methylated CpG and
# the occurences of that profile in the sample
# INPUT:
# sfiles   -> list of .out file paths containing methylation profiles of each sample generated with getProfilesFromBlast.py module
# out_name -> output file name
def makeAbTab(sfiles, out_name):
    sample = []
    sample += [re.sub(".out", "", ntpath.basename(name).split('_', 1)[-1]) for name in sfiles]

    nSam = len(sample)

    tab_out = open(out_name, "w")

    profiles = {}
    for i in range(nSam):

        file = open(sfiles[i])
        # counts occurences of each methylation profile in the sample
        for line in file:
            line = line.rstrip()
            if not line.__contains__("2"):
                if line in profiles:
                    profiles[line][i] += 1
                else:
                    profiles[line] = [0] * nSam
                    profiles[line][i] = 1

    keys = sorted(profiles.keys())
    nSpec = len(keys)

    tab_out.write('profile\tn_meths\t' + "\t".join(sample) + '\n')

    for i in xrange(nSpec):
        n_met = sum([int(val) for val in keys[i].replace(" ", "")])
        tab_out.write(keys[i].replace(" ", "") + '\t' + str(n_met) +
                      '\t' + "\t".join([str(n) for n in profiles[keys[i]]]) + "\n")

    tab_out.close()


# for each analyzed region builds summary files containing profile aboundances
# two files per region are produced: one in tabular form and another in BIOM format
# INPUT:
# samples -> list of samples IDs
# regions -> list of regions IDs
# outDir  -> output folder path
def make_summary(samples, regions, outDir):
    samPaths = []
    for reg in regions:
        destDir = outDir + "/" + reg
        os.mkdir(destDir)
        for sam in samples:
            samPaths.append(sam + "/" + reg + "_" + sam + "/" + reg + "_" + sam + ".out")
        makeAbTab(samPaths, destDir + "/" + reg + ".csv")
        makeBiom(samPaths, destDir + "/" + reg + ".biom")
        samPaths = []


# performs a series of qualitative analises using local qiime installation for each one of the analyzed regions
# parameters
# INPUT:
#   biom_reg        -> biom file containing methylation profile abundances for each sample
#   destDir         -> output folder
#   reg             -> id of the analysed region
#   sampleMeta      -> tab separated file containing samples to analyse along with sample features
#   qiime_env_path  -> path to local qiime environment script
def qualitative_Analysis(biom_reg, destDir, reg, sampleMeta, qiime_env_path):
    # count number of samples to process from the metadata file
    nSam = 0
    f = open(sampleMeta)
    for line in f:
        nSam += 1
    nSam -= 1
    # sort samples based on they Description field
    os.system(".  " + qiime_env_path + " && " + "sort_otu_table.py -i " + biom_reg + " -m" + sampleMeta +
              " -s Description -o " + destDir + "/" + reg + "_sorted.biom >> " + destDir + "/qiime.log 2>&1")
    # compute distribution of the total number of methylation profiles for the given set of samples
    os.system(
        "biom summarize-table -i " + biom_reg + " -o " + destDir + "/" + reg + "_summary.txt >> " + destDir + "/qiime.log 2>&1 ")


    # prepare alfa and beta params files for quiime analyses

    f_out = open("alpha_params.txt", "w")
    f_out.write("alpha_diversity:metrics observed_species,shannon,chao1,simpson,singles"+"\n"+"summarize_taxa:level 1"+"\n"+"plot_taxa_summary:chart_type bar,pie")
    f_out.close()

    f_out = open("beta_params.txt", "w")
    f_out.write("beta_diversity:metrics  bray_curtis")
    f_out.close()


    # summarize samples composition based on the number of methylated sites per profile
    os.system(
        ".  " + qiime_env_path + " && " + "summarize_taxa_through_plots.py -f -i " + biom_reg + " -m" + sampleMeta +
        " -p alpha_params.txt -o " + destDir + "/profileSummary/ >> " + destDir + "/qiime.log 2>&1 ")
    f = open(destDir + "/" + reg + "_summary.txt")
    l = f.readline()
    while l and not l.startswith(" Min"):
        l = f.readline()
    depth = str(int(float(l.rstrip().split(":")[1])))
    # if there is more than one sample perform comparative analyses among samples (Alpha and Beta diversities)
    if nSam > 1:
        # build an heatmap with profile aboundaces in each sample
        os.system(". " + qiime_env_path + " && " + "make_otu_heatmap.py -i " + biom_reg + " -m" + sampleMeta +
                  " -o " + destDir + "/heatmap.pdf >> " + destDir + "/qiime.log 2>&1")
        # compute Alpha diversity metrics reported in alpha_params.txt file with a rarefaction procedure
        os.system(". " + qiime_env_path + " && " + "alpha_rarefaction.py -i " + biom_reg + " -m" + sampleMeta +
                  " -p alpha_params.txt -o " + destDir + "/Alpha>> " + destDir + "/qiime.log 2>&1")
        # compute Beta diversity analyses based on PCA and bray-curtis distances among samples
        os.system(
            ". " + qiime_env_path + " && " + "beta_diversity_through_plots.py -i " + biom_reg + " -m" + sampleMeta +
            " -p beta_params.txt -o " + destDir + "/Beta -e " + depth + " >> " + destDir + "/qiime.log 2>&1")
        # make 2D plots of samples witgh principal components from the previous analisys
        os.system(
            ". " + qiime_env_path + " && " + "make_2d_plots.py -i " + destDir + "/Beta/bray_curtis_pc.txt -m" +
            sampleMeta + " -o " + destDir + "/Beta/PCA >> " + destDir + "/qiime.log 2>&1")
        # make distance boxplots based on bray-curtis distance among samples
        os.system(
            ". " + qiime_env_path + " && " + "make_distance_boxplots.py -d " + destDir + "/Beta/bray_curtis_dm.txt -m" +
            sampleMeta + " -f Description -o " + destDir + "/Beta/dist_boxplot >> " + destDir + "/qiime.log 2>&1")
