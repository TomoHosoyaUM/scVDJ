#!/bin/bash
# usage:   bash 21_fa_ROW_COL_DV.sh

DATEstart=$(date),  UTCstart=$(date +%s)
echo "### $0 ###"                                                          >> ./analysis.log.txt
echo "### $DATEstart ------------------------------ ### START ### $0 ###"  >> ./analysis.log.txt
echo "### $DATEstart ------------------------------ ### START ### $0 ###"  

#################################################################################################
## input fasta ##################################################################################
input=../input_CCS/*.fa

find ${input}                                           >> ./analysis.log.txt
find ${input}

Nseq=$(cat ${input} | grep ">" | wc -l )
echo "### ${Nseq} total CCSs in the input FASTA files"  >> ./analysis.log.txt
echo "### ${Nseq} total CCSs in the input FASTA files"  

## starting with barcode-AF ####################################################################
echo "### fa_ROW ### starting with barcode-AF"
mkdir -p ./fa_ROW

bcfile="../barcode/bc_row.txt"
i="1"
while read -r line
do
    barcode="${line}"
    cat ${input} | grep -B1 "${barcode}" | seqkit seq -w 0 -rp > ./fa_ROW/ROW${i}.rc.fa  # save reverse-complement for now
let "i++"
done < "${bcfile}"

## ending with barcode-AF ######################################################################
echo "### fa_ROW ### ending with barcode-AF"
bcfile="../barcode/bc_row_rc.txt"
i="1"
while read -r line
do
    barcode="${line}"
    cat ${input} | grep -B1 "${barcode}" | sed -e 's/>/>rc./g' >> ./fa_ROW/ROW${i}.rc.fa
    
    # remove after barcode-AF
    # from .........barcode-AF(rc)..........
    # to   .........barcode-AF(rc)
    # If two barcode-AF(rc) in a sequence, only first one will be kept
    cat fa_ROW/ROW${i}.rc.fa | sed -e "s/${barcode}.*/${barcode}/g" | seqkit seq -w 0 -rp > ./fa_ROW/ROW${i}.fa
let "i++"
done < "${bcfile}"

# remove ./fa_ROW/*.rc.fa
rm ./fa_ROW/*.rc.fa

## ending with barcode-AR3 in "./fa_ROW/ROW*.fa #################################################
echo "### fa_ROW_COL"
mkdir -p ./fa_ROW_COL

bcfile="../barcode/bc_col_rc.txt"
i="1"
while read -r line
do
    barcode="${line}"
    echo "${barcode}"
    for FILE in ./fa_ROW/ROW*.fa
        do
        cat ${FILE} | grep -B1 "${barcode}" | sed -e "s/${barcode}.*/${barcode}/g" | sed -e 's/-//g' > ./fa_ROW_COL/$(basename ${FILE} .fa).COL${i}.fa
    done
let "i++"
done < "${bcfile}"

# remove directory ./fa_ROW
rm -r -f ./fa_ROW

## number of CCS in each well (ROW.COL) #########################################################
echo "### calculate number of CCSs in each well"
OUT=./table21_ROW_COL_full.csv
echo -n > ${OUT}

for i in `seq 1 200`
do
    for j in `seq 1 12`
    do
        barcode=ROW${i}.COL${j}
        NUM=$(cat ./fa_ROW_COL/ROW${i}.COL${j}.fa | grep '>' | wc -l)
        echo "${barcode},${NUM}" >> ${OUT}
    done
done
cat $OUT | awk 'BEGIN{OFS=","} $2 >0' > ./table21_ROW_COL_short.csv

## V or D primers ###############################################################################
echo "### fa_ROW_COL_DV"
mkdir -p fa_ROW_COL_DV

bcfile="../barcode/DV_seq.txt"
i="1"
while read -r line
do
    barcode="${line}"
    echo "${barcode}"
    for FILE in fa_ROW_COL/*.fa
    do
        # barcode 14 bp, AF 18 bp, DV, 17-25 bp
        seqkit grep -w 0 -s -R 33:57 -r -i -p "${barcode}" ${FILE} > ./fa_ROW_COL_DV/$(basename $FILE .fa).DV${i}.fa
    done
let "i++"
done < "${bcfile}"

# remove empty
find ./fa_ROW_COL_DV -type f -empty -delete

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
