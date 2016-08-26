#!/usr/local/bin/python

################################################################################
#
# check_length.py: checks the length of input multi-FASTA file. Reports
#                  the average length, sd, and range
#
# Required Inputs: [1] input file name
#
# written by Mario Muscarella
# last update 26 Aug 2016
#
################################################################################

import sys
import os
from Bio import SeqIO
import numpy

in_file = sys.argv[1]

seq_len = list()

for seq_record in SeqIO.parse(in_file, "fasta"):

    seq_len_raw = len(seq_record.seq)
    # print(seq_len_raw)
    seq_len.append(seq_len_raw)

sum(seq_len)/len(seq_len)
avg = numpy.mean(seq_len)
sd = numpy.std(seq_len)
mins = numpy.amin(seq_len)
maxs = numpy.amax(seq_len)

print("average sequence length: %i" % (avg))
print("standard deviation: %i" % (sd))
print("range: [%i, %i]" % (mins, maxs))
