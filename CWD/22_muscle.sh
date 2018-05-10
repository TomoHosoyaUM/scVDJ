#!/bin/bash
# usage:   bash 22_muscle.sh

DATEstart=$(date),  UTCstart=$(date +%s)
echo "### $0 ###"                                                          >> ./analysis.log.txt
echo "### $DATEstart ------------------------------ ### START ### $0 ###"  >> ./analysis.log.txt
echo "### $DATEstart ------------------------------ ### START ### $0 ###"  

#################################################################################################
# muscle alignment
mkdir -p ./muscle_fastaout
mkdir -p ./muscle_htmlout
mkdir -p ./muscle_fastaout_nowrap

for FILE in ./fa_ROW_COL_DV/ROW*.fa
do
    muscle -in $FILE -fastaout ./muscle_fastaout/$(basename $FILE .fa).afa -htmlout ./muscle_htmlout/$(basename $FILE .fa).html
done

echo "nowrap.afa"
for FILE in ./muscle_fastaout/ROW*.afa
do
    cat $FILE | seqkit seq -w 0 > ./muscle_fastaout_nowrap/$(basename $FILE .afa).fa
done

#################################################################################################
# BASH DONE
DATEend=$(date),  UTCend=$(date +%s)
i=$(echo "${UTCend} - ${UTCstart}" | bc )
((sec=i%60, i/=60, min=i%60, hrs=i/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "### ${DATEstart} - ${DATEend} ### DONE ### $0 ###"
echo "### ${DATEstart} - ${DATEend} ### DONE ### $0 ###"  >> ./analysis.log.txt
echo "### ${timestamp}   hh:mm:ss"
echo "### ${timestamp}   hh:mm:ss"                        >> ./analysis.log.txt
echo " "                                                  >> ./analysis.log.txt
