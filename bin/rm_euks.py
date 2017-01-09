import sys
import os
import re
from Bio import SeqIO
import numpy
import pandas

#in_file = sys.argv[1]
#out_file = sys.argv[2]
in_file = "./data/nmicrobiol201648-s7.txt"
out_file = "./data/NTL.clean.fasta"

NonEuks = []
names_all = []
names_noneuks = []

for seq_record in SeqIO.parse(in_file, "fasta"):
    temp = seq_record.name
    names_all.append(temp)
    if not re.match("Euk", temp):
        NonEuks.append(seq_record)
        names_noneuks.append(temp)

SeqIO.write(NonEuks, out_file, "fasta")

num_euks = len(names_all) - len(names_noneuks)
num_noneuks = len(names_noneuks)

print("Done")
print("Removed %i Eukaryote sequences and saved %i Bacteria and Archaea") % (num_euks, num_noneuks)
