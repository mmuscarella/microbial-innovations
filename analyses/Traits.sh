#!/bin/bash

rm ../data/LTPs123.faprotax
rm ../data/report.txt

../bin/make_table.py ../data/LTPs123_unique.nr_v123.wang.taxonomy ../data/LTPs123.OTUtable
../bin/collapse_table.py -i ../data/LTPs123.OTUtable -g ../bin/FAPROTAX.txt -o ../data/LTPs123.faprotax --omit_columns 0 -d 'taxonomy' -c '#' -v -r ../data/report.txt
