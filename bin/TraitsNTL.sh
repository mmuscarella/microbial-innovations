#!/bin/bash

rm ../data/NTL.faprotax
rm ../data/NTLreport.txt

../bin/make_table.py ../data/NTL.clean.pds.wang.taxonomy ../data/NTL.OTUtable
../bin/collapse_table.py -i ../data/NTL.OTUtable -g ../bin/FAPROTAX.txt -o ../data/NTL.faprotax --omit_columns 0 -d 'taxonomy' -c '#' -v -r ../data/NTLreport.txt

rm ../data/NTL.OTUtable
