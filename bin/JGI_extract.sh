#!/bin/bash

################################################################################
#
# JGI_extract.sh: Runs
#
# Required Inputs:
#
# written by ME Muscarella
# last update 2017 04 17
#
################################################################################


# Extract Archive Files
cd ~/JGI/output

files=$( ls ./img/*.tar.gz )
echo "$files"

for i in $files ; do
  echo $i
  name=$( echo $i | grep -o "\.[0-9]*\." | sed $'s/\.//g' )
  if [[ ! -d $name ]]; then
    tar -xkzvf $(echo $i )
  fi
done

# Extract KEGG KO
Rscript ~/microbial-innovations/bin/IMG_KO.R

# Extract FASTA Sequences
gff=$( ls ./*/*.gff )
total_partial=0
total_full=0
cat /dev/null > JGI.fasta
echo -e "Genome\tFullLength\tPartialLength" > JGIseqs.txt

for i in $gff ; do
  echo $i
  genome=$( echo $i | grep -o "/.*/" | sed $'s|\/||g' )
  partial=0
  full=0
  rRNA=$( grep "rRNA.*product=16S" $(echo $i) )
  rRNA=$( echo "$rRNA" | sed 's/^.*ID=//' | sed 's/\;locus.*//' )
  echo "$rRNA"
  fna=$( echo $i | sed 's/.gff/.genes.fna/' )
  if [ ! -e $fna ] ; then
    echo "Not Found"
    fna=$( ls $( echo $i | sed 's/.gff/\*.fna/' ) )
  fi 
  echo "$fna"
  for j in $rRNA ; do
    id=$( grep $( echo $j ) $fna )
    echo "$id"
    seq=$( sed -e '/'"$j"'/,/^\s*$/!d' $fna )
    echo "$seq"
    echo "$seq" | wc -c
    echo "$id" | wc -c
    echo "$seq" | sed -n '1!p' | wc -c
    seq_len=$( echo "$seq" | sed -n '1!p' | wc -c )
    if [[ $seq_len > 1200 ]] ; then
      echo Yes
      total_full=$(($total_full + 1))
      full=$(($full + 1))
      echo "$seq" >> JGI.fasta
    else
      echo NO
      total_partial=$((total_partial + 1))
      partial=$((artial + 1))
    fi
  done
  echo "Found $full full length sequences"
  echo "Found $partial partial length sequences"
  echo -e "$genome\t$full\t$partial" >> JGIseqs.txt

done
echo "16S rRNA gene search complete"
echo "Found $total_full full length sequences"
echo "Found $total_partial partial length sequences"
