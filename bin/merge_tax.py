import sys
import os
import re
from Bio import SeqIO
import numpy
import pandas

in_file1 = sys.argv[1]
in_file2 = sys.argv[2]
out_file = sys.argv[3]
# in_file1 = './data/LTPs123_unique.nr_v123.wang.taxonomy'
# in_file2 = './data/LTPs123_unique.taxonomy'
# out_file = './data/LTPs123.unique.full.taxonomy'

tax1 = pandas.DataFrame.from_csv(in_file1, sep = '\t', header = None, index_col = None)
tax2 = pandas.DataFrame.from_csv(in_file2, sep = '\t', header = 0, index_col = None)

if(len(tax1) == len(tax2)):
    len_tax = len(tax1)
else:
    print('Stop: Taxonomy files are not the same length')
    sys.exit(1)

if not(numpy.array_equal(tax1.iloc[0:, 0], tax2.iloc[0:, 0])):
    print('Stop: OTUs are not ordered correctly')
    sys.exit(1)

print(tax1.iloc[0:3, 0:])
print(tax2.iloc[0:3, 0:])

temp_otu = list()
temp_tax = list()

for i in range(0,len_tax):
    otu = tax1.iloc[i, 0]
    temp_otu.append(otu)
    temp = tax1.iloc[i,1]
    temp_split = temp.rsplit(';')
    temp_high = '%s;%s;%s;%s' % (temp_split[0], temp_split[1], temp_split[2], temp_split[3])
    temp_low = tax2.iloc[i, 1]
    temp_join = '%s;%s' % (temp_high, temp_low)
    temp_tax.append(temp_join)

OTU_Tax = pandas.DataFrame(
    {'OTU': temp_otu,
     'Taxonomy': temp_tax
    }, columns = ['OTU', 'Taxonomy'])

OTU_Tax.to_csv(out_file, sep='\t', header=False, index=False)

print("Done")
