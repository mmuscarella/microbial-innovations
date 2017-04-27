#!/bin/bash

################################################################################
#
# JGI.sh: Runs
#
# Required Inputs:
#
# written by ME Muscarella
# last update 2017 04 12
#
################################################################################

sed $'s/\t\t/\tNA\t/g' taxontable64227_17-apr-2017.xls | sed $'s/\t\t/\tNA\t/g' | sed $'s/\t\t/\tNA\t/g' > bacteriataxontable.txt

sed $'s/\t\t/\tNA\t/g' taxontable32958_18-apr-2017.xls | sed $'s/\t\t/\tNA\t/g' | sed $'s/\t\t/\tNA\t/g' > archaeataxontable.txt

cat bacteriataxontable.txt > taxontable.txt
awk FNR-1 archaeataxontable.txt >> taxontable.txt

out="tee -a JGI_DB.logfile"

while IFS=$'\t' read -r -a myArray ; do
  JGI=${myArray[16]}
  if [[ $JGI == 'JGI Project ID / ITS PID' ]] ; then
    echo "The Mission to Conquer The Microbial World Has Begun" | tee JGI_DB.logfile
    echo -e "${myArray[4]}\t${myArray[16]}\tDB\tRel_DB" > JGImatch.txt
  else
    curl -i 'http://genome.jgi.doe.gov/'$JGI'/status' > temp1.txt
    cat temp1.txt
    db=$(grep -o 'db=.*' temp1.txt | sed 's/db=//')
    db=${db%$'\r'}
    echo "Genome analysis ID: $db" | $out

    curl "http://genome.jgi.doe.gov/$db/$db.info.html" > temp2.txt
    rel_url=$(grep -o "jgiProjectId.*Info" temp2.txt)
    rel=$( echo $rel_url | grep -o "keyValue=.*&amp;"  | sed $'s/keyValue=//' | sed $'s/\&amp\;//' )
    echo "Related Project ID: $rel" | $out

    curl -i 'http://genome.jgi.doe.gov/'$rel'/status' > temp3.txt
    cat temp3.txt
    db2=$(grep -o 'db=.*' temp3.txt | sed 's/db=//')
    if [[ ${db2%$'\r'} == "db=null" ]] ; then
      db2=NA
    else
      echo "Related analysis ID: $db2" | $out
      db2=${db2%$'\r'}
    fi

    echo -e "${myArray[4]}\t${myArray[16]}\t$db\t$db2" >> JGImatch.txt
    rm ./temp1.txt
    rm ./temp2.txt
    rm ./temp3.txt
  fi

  # Sleep to prevent overdoing it on server
  sleep 10

done < ./taxontable.txt

# Print Result to Screen
result=$( awk 'BEGIN { FS = "\t" } ; { print $3 }' JGImatch.txt )
result2=$( awk 'BEGIN { FS = "\t" } ; { print $4 }' JGImatch.txt )
total=$( echo "$result" | wc -l )
colname=$( echo "$result" | grep "DB" | wc -l )
fail=$( echo "$result" | grep "null" | wc -l )
fail2=$( echo "$result2" | grep "null" | wc -l )
echo " $(($total - $fail - $colname)) out of $(($total - $colname)) Genomes projects were found"
echo " $(($total - $fail2 - $colname)) out of $(($total - $colname)) Genomes are part of larger projets"
