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
cat /dev/null > JGI.first.fasta
cat /dev/null > JGI.rename.fasta
cat /dev/null > JGI.full.rename.fast[a
cat /dev/null > JGI.first.rename.fasta
echo -e "Genome\tFullLength\tPartialLength" > JGIseqs.txt
echo -e "Name\tGenomeName\tGenome" > JGIseqMatch.txt
echo -e "Name\tGenomeName\tGenome" > JGIseqMatchfirst.txt

for i in $gff ; do
  echo $i
  genome=$( echo $i | grep -o "/.*/" | sed $'s|\/||g' )
  partial=0
  full=0
  rRNA=$( grep "rRNA.*product=16S" $(echo $i) )
  rRNA=$( echo "$rRNA" | sed 's/^.*ID=//' | sed 's/\;locus.*//' | sed 's/\;[a-zA-Z].*//' )
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
    if [[ -z $id ]] ; then
      echo "Not Found"
      j=$( echo $j | sed 's/\.[0-9]*//' )
      id=$( grep $( echo $j ) $fna )
      echo "$id"
    fi
    name=$( echo $iqd | cut -d " " -f 2 )
    echo "$name"
    seq_entry=$( sed -e '/'"$j"'/,/^\s*$/!d' $fna )
    if (( $( echo $seq_entry | wc -c ) > 5000 )); then
      echo "Possibly Full Genome; Try ending with >"
      seq_entry=$( sed -e '/'"$j"'/,/>/!d' $fna | sed '$d' )
      name=$( echo $id | sed 's/>//' )
    fi
    echo "$seq_entry"
    seq=$( echo "$seq_entry" | sed -n '1!p' )
    echo "$seq" | wc -c
    echo "$id" | wc -c
    seq_len=$( echo "$seq" | wc -c )
    new_name=">$j $genome $name"
    echo "$new_name" >> JGI.full.rename.fasta
    echo "$seq" >> JGI.full.rename.fasta
    if (( $seq_len > 1199 )) ; then
      echo Yes
      total_full=$(($total_full + 1))
      full=$(($full + 1))
      echo "$new_name" >> JGI.fasta
      echo "$seq" >> JGI.fasta
      echo ">$genome" >> JGI.rename.fasta
      echo "$seq" >> JGI.rename.fasta
      echo -e "$j\t$name\t$genome" >> JGIseqMatch.txt
      if (( $seq_len > 1199 && $full == 1 && $seq_len < 2000)); then
        echo Yes
        echo "$new_name" >> JGI.first.fasta
        echo "$seq" >> JGI.first.fasta
        echo ">$genome" >> JGI.first.rename.fasta
        echo "$seq" >> JGI.first.rename.fasta
        echo -e "$j\t$name\t$genome" >> JGIseqMatchfirst.txt
      fi
    else
      echo "No Full Match"
      total_partial=$((total_partial + 1))
      partial=$((partial + 1))
    fi
  done
  echo "Found $full full length sequences"
  echo "Found $partial partial length sequences"
  echo -e "$genome\t$full\t$partial" >> JGIseqs.txt

done
echo "16S rRNA gene search complete"
echo "Found $total_full full length sequences"
echo "Found $total_partial partial length sequences"
