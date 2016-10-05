#!/usr/local/bin/python

################################################################################
#
# name_change.py: renames RDP database to replace the sequence ID with the
#                 appropriate accession number which can be used for NCBI
#                 queries
#
# Required Inputs: [1] input file name [2] output file name
#
# written by Mario Muscarella
# last update 30 Sep 2016
#
################################################################################

import sys
import os
import re
from Bio import SeqIO

# in_file = "./"
# out_file ="./"

in_file = sys.argv[1]
out_file = sys.argv[2]

output_handle = open(out_file, "w")

OldseqID = list()
seqInfo = list()
NewseqID = list()

for seq_record in SeqIO.parse(in_file, "fasta"):
    OldseqID.append(seq_record.id)
    seqInfo.append(seq_record.description)
    temp = seq_record.description
    temp_split = temp.rsplit('; ', 1)
    temp_final = temp_split[1]
    seqID_New = temp_final
    NewseqID.append(seqID_New)
    seq_record.id = seqID_New

    SeqIO.write(seq_record, output_handle, "fasta")

output_handle.close()
renamed = len(NewseqID)

print("done")
print("renamed %i sequence records" %(renamed))
