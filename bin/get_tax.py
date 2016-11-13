import sys
import os
import re
from Bio import SeqIO
import numpy
import pandas

in_file = sys.argv[1]
out_file = sys.argv[2]
# in_file = "./data/LTPs123_unique.fasta"
# out_file = "./data/LTPs123.unique.taxonomy"

OTU = list()
Tax1 = list()
Tax2 = list()
Tax3 = list()
Tax4 = list()

for seq_record in SeqIO.parse(in_file, "fasta"):
    OTU.append(seq_record.id)
    temp = seq_record.description
    temp_split = temp.rsplit("\t")
    Tax1.append(temp_split[6])
    temp2 = temp_split[5]
    temp2_split = temp2.rsplit(" ")
    Tax2.append(temp2_split[0])
    Tax3.append(temp2_split[1])
    temp3 = "%s;%s;%s" % (temp_split[6], temp2_split[0], temp2_split[1])
    Tax4.append(temp3)

OTU_Tax = pandas.DataFrame(
    {'OTU': OTU,
     'Family': Tax1,
     'Genus': Tax2,
     'Species': Tax3
    }, columns = ['OTU', 'Family', 'Genus', 'Species'])

OTU_Tax2 = pandas.DataFrame(
    {'OTU': OTU,
     'Taxonomy': Tax4
     }, columns = ['OTU', 'Taxonomy'])

OTU_Tax2.to_csv(out_file, sep='\t', header=True, index=False)

print("done")
