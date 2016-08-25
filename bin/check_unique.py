#!/usr/local/bin/python

################################################################################
#
# check_unique.py: checks an input multifasta to determine if record IDs are
#                  uniuqe. If IDs are not unique the script uses only the first
#                  instance of a given ID and saves a new "uniuqe" output.
#
# Required Inputs: [1] input file name [2] output filen name
#
# written by Mario Muscarella
# last update 25 Aug 2016
#
################################################################################

import sys
import os
from Bio import SeqIO

if len(sys.argv) > 1:
    in_file = sys.argv[1]

    if len(sys.argv) > 2:
        out_file = sys.argv[2]
    else:
        out_file_raw = os.path.splitext(in_file)[1]
        out_file = out_file_raw+".unique.fasta"
else:
    print("you must supply at least an input file")

seqs = list()

for seq_record in SeqIO.parse(in_file, "fasta"):

    seqs.append(seq_record.id)

diff_id = list({x for x in seqs if seqs.count(x) > 1})
diff_len = len(diff_id)
diff_ids = str(diff_id)[1 : -1]

print("There are %i replicate IDs: %s" % (diff_len, diff_ids))

if diff_len > 0:

    output_handle = open(out_file, "w")

    seq_ids = list()
    seq_saved = list()

    for seq_record in SeqIO.parse(in_file, "fasta"):
        seq_id = seq_record.id
        seq_ids.append(seq_id)

        if seq_id not in diff_id:
            seq_saved.append(seq_id)
            count = SeqIO.write(seq_record, output_handle, "fasta")

    output_handle.close()

    seq_orig = len(seq_ids)
    seq_save = len(seq_saved)

    print("Saved %i of %i FASTA records" % (seq_save, seq_orig))
