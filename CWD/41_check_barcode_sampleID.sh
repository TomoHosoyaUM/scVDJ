#!/bin/bash
echo "### $0 ###"                        >> ./analysis.log.txt
echo "### $0 ###" 
# usage:   bash 41_check_input_files.sh

# check ../barcode_sampleID/IDcell_*.txt
echo "### Number of barcode ROW.COL"     >> ./analysis.log.txt
echo "### Number of barcode ROW.COL"
wc -l ../barcode_sampleID/IDcell_*.txt   >> ./analysis.log.txt
wc -l ../barcode_sampleID/IDcell_*.txt 

cat ../barcode_sampleID/IDcell_*.txt | sort | uniq -d > ../barcode_sampleID/barcode_duplicate.txt
Ndup=$(cat ../barcode_sampleID/barcode_duplicate.txt | wc -l )

if [ ${Ndup} -gt 0 ]; then 
    echo "### ${Ndup} duplicate found for barcodes in  ../barcode_sampleID/IDcell_*.txt" >> ./analysis.log.txt
    echo "### ${Ndup} duplicate found for barcodes in  ../barcode_sampleID/IDcell_*.txt"
    echo "### ${Ndup} duplicate list is            in  ../barcode_sampleID/barcode_duplicate.txt"
fi

# BASH DONE
