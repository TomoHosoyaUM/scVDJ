#!/bin/bash
# usage:   bash 46_sampleID.sh

echo "### $0 ###"   >> ./analysis.log.txt
echo "### $0 ###" 

##########################################################################################
OUTdir="./table45_sampleID"
mkdir -p ${OUTdir}

Rscript 45a_remove_error.r

cd ${OUTdir}
Rscript ../45b_sampleID.r ../../barcode_sampleID/IDcell_*.txt
cd ..

##########################################################################################
# BASH DONE
