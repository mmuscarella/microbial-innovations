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

read -s -p "Please Enter Password: `echo $'\n> '`" my_password

curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=mmuscar@illinois.edu' --data-urlencode "password=$my_password" -c cookies > /dev/null

out="tee -a JGI_Download.logfile"

if [ ! -d ~/JGI ]; then
  mkdir ~/JGI
fi

cd ~/JGI

if [ ! -d output ]; then
  mkdir output
fi

if [ ! -d output/gbk ]; then
  mkdir output/gbk
fi

if [ ! -d output/img ]; then
  mkdir output/img
fi

while IFS=$'\t' read -r -a myArray ; do

  JGI=${myArray[1]}

  if [[ $JGI == 'JGI Project ID / ITS PID' ]] ; then
    echo "The Mission to Conquer The Microbial World Has Begun"  | tee JGI_Download.logfile
    echo -e "${myArray[0]}\t${myArray[1]}\t${myArray[2]}\tIMG_Name\tIMG_Downloaded" > JGIdownload.txt
  else
    DB1=${myArray[2]}
    DB1=${DB1%$'\r'}
    DB2=${myArray[3]}
    DB2=${DB2%$'\r'}
    echo $JGI
    echo $DB1
    echo $DB2

    db1=$(curl -i 'http://genome.jgi.doe.gov/'$DB1'/status' | grep -o "db=.*")
    db2=$(curl -i 'http://genome.jgi.doe.gov/'$DB2'/status' | grep -o "db=.*")

    if [[ ${db2%$'\r'} == "db=null" ]] ; then
      curl "http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism=$DB1" -b cookies > temp.xml
      DB=$DB1
      else
      curl "http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism=$DB2" -b cookies > temp.xml
      DB=$DB2
    fi

    sed -e $'/folder name\=\"IMG Data\"/,/folder/!d' temp.xml > temp2.xml
    #gbk=$(grep '.gbk' temp2.xml)
    data=$(grep 'tar.gz' temp2.xml)

    #if [ ! -z "$gbk" ] ; then
    #  echo "$gbk" | $out
    #  gbk_result=YES
    #
    #  if [[ $( echo "$gbk" | wc -l ) == 1 ]] ; then
    #    echo "Single GBK Entry" | $out
    #    gbk_url=$( echo "$gbk" | grep -o 'url.*gbk\"' | sed $'s/url=\"//' \
    #                        | sed $'s/\"//' | sed $'s/\&amp;/\&/' )
    #    gbk_url=${gbk_url%$'\r'}
    #    gbk_name="$DB.$( echo $gbk_url | rev | cut -d'/' -f 1 | rev )"
    #    curl "http://genome.jgi.doe.gov$gbk_url" -b cookies > ./output/gbk/$gbk_name
    #  else
    #    echo "Multiple GBK Entries" | $out
    #    #timestamp=$( echo "$gbk" | grep -o -P 'timestamp=\"(.*?)\"' )
    #    gbk_url=$( echo "$gbk" | grep -o 'url.*gbk\"' | sed $'s/url=\"//' \
    #                        | sed $'s/\"//' | sed $'s/\&amp;/\&/' )
    #    gbk_url=$( echo "$gbk_url" | head -n1 )
    #    gbk_url=${gbk_url%$'\r'}
    #    gbk_name="$DB.$( echo $gbk_url | rev | cut -d'/' -f 1 | rev )"
    #    curl "http://genome.jgi.doe.gov$gbk_url" -b cookies > ./output/gbk/$gbk_name
    #  fi
    #else
    #  echo "No Genbank File Found for $DB"  | $out
    #  gbk_result=NO
    #fi

    if [ ! -z "$data" ] ; then
      echo $data
      data_result=YES

      if [[ $( echo "$data" | wc -l ) == 1 ]] ; then
        echo "Single IMG Entry" | $out
        data_url=$( echo $data | grep -o 'url.*tar.gz\"' | sed $'s/url=\"//' \
                          | sed $'s/\"//g' | sed $'s/\&amp;/\&/' )
        data_url=${data_url%$'\r'}
        data_name="$DB.$( echo $data_url | rev | cut -d'/' -f 1 | rev )"
        if [[ ! -e ./output/img/$data_name ]] ; then
          echo "Not Found"
          curl "http://genome.jgi.doe.gov$data_url" -b cookies > ./output/img/$data_name
        fi

      else
        echo "Multiple IMG Entries" | $out
        #timestamp=$( echo "$gbk" | grep -o -P 'timestamp=\"(.*?)\"' )
        data_url=$( echo $data | grep -o 'url.*tar.gz\"' | sed $'s/url=\"//' \
                          | sed $'s/\"//g' | sed $'s/\&amp;/\&/' )
        data_url=$( echo "$data_url" | head -n1 )
        data_url=${data_url%$'\r'}
        data_name="$DB.$( echo $data_url | rev | cut -d'/' -f 1 | rev )"
        if [ ! -e ./output/img/$data_name ] ; then
          echo "Not Found"
          curl "http://genome.jgi.doe.gov$data_url" -b cookies > ./output/img/$data_name
        fi
      fi
    else
      echo "No IMG Files Found for $DB"  | $out
      data_name=NA
      data_result=NO
    fi
    rm ./temp.xml
    rm ./temp2.xml
    echo -e "${myArray[0]}\t${myArray[1]}\t$DB\t$data_name\t$data_result" >> JGIdownload.txt
  fi

  # Sleep to prevent overdoing it on server
  sleep 10

done < ./JGImatch.txt

# Print Result to Screen
result=$( awk 'BEGIN { FS = "\t" } ; { print $5 }' JGIdownload.txt )
total=$( echo "$result" | wc -l )
colname=$( echo "$result" | grep "IMG_Downloaded" | wc -l )
complete=$( echo "$result" | grep "YES" | wc -l )
fail=$( echo "$result" | grep "NO" | wc -l )
echo " $complete out of $(($total - $colname)) Genomes Downloaded"
