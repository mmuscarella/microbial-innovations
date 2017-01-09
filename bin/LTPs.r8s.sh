#!/bin/bash

rm ../data/LTP.r8s.txt
cat ../data/LTP.rooted.nex ../bin/r8s.config > ../data/LTP.r8s.txt

r8s -b -f ../data/LTP.r8s.txt
