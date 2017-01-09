#!/bin/bash

echo "Parsing Green Genes Data"

rm ./IMG.fasta

cat GG_to_IMGv350.txt | while read -a p; do
  IMG=${p[1]}
  echo $IMG
  IMGGG=$(grep "$IMG" gg_13_5_img.txt)
  echo $IMGGG
  GG=( $IMGGG )
  #gg=${GG[0]}
  #echo $gg
  #echo $gg >> test.txt
  grep -m1 -A1 "^>$GG$" gg_13_5.fasta >> IMG.fasta
  # grep -A1 "^>$gg$" gg_13_5.fasta | wc -l
  # that part could be improved, search for best sequence
  # sometimes it didn't find a sequence at all and that sucks
  # also, it is a bit slow
  len=$(grep '>' IMG.fasta | wc -l)
  echo $len
done
