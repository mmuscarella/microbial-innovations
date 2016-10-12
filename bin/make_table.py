#!/usr/local/bin/python

################################################################################
#
# make_table.py: Takes a list of bacterial strains or OTUs and makes an OTU
#                table where each site has only one strain with abundance one.
#                The purpose of this is to prepare a bacterial list for
#                function calling programs such as FAPROTAX.
#
# Required Inputs: [1] input file name as a taxonomy file
#                  [2] output file name
#
#                  Taxonomy file should be mothur output based on SILVA
#
# written by Mario Muscarella
# last update 11 Oct 2016
#
################################################################################

import sys
import os
from Bio import SeqIO
import numpy
import pandas

in_file = sys.argv[1]
out_file = sys.argv[2]

# Test
# in_file = './data/LTPs123_unique.pds.wang.taxonomy'
# out_file = './data/test.OTUtable'


if in_file.endswith('.taxonomy'):
    tax = pandas.DataFrame.from_csv(in_file, sep='\t', header = None, index_col=False)
else:
    print('Input must be a Taxonomy file')

name = tax.loc[0:, 0].values
taxonomy = tax.loc[0:, 1].values
par = len(name)
print ("Table includes %s OTUs" % par)

taxa = pandas.DataFrame(numpy.zeros((par, par)), columns = name, index = name).astype(int)

if (taxa.columns.values[0:10] == taxa.index.values[0:10]).all():
    print("success")
    taxa.values[[numpy.arange(par)]*2] = 1
else:
    for i in name:
        print i
        for j in name:
            if i == j:
                taxa.set_value(i,j,1)
            else:
                taxa.set_value(i,j,0)

taxa.insert(0, 'taxonomy', taxonomy)

# print taxa.iloc[0:5,0:5]
taxa.to_csv(out_file, sep='\t', header = True, index = True)
print("done")
