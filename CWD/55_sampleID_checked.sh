#!/bin/bash
# usage:   bash 55_sampleID_checked.sh

echo "### $0 ###"   >> ./analysis.log.txt
echo "### $0 ###" 

##########################################################################################
OUTdir="./table55_sampleID_checked"
mkdir -p ${OUTdir}

Rscript 55a_remove_checked.r

cd ${OUTdir}
Rscript ../45b_sampleID.r ../../barcode_sampleID/IDcell_*.txt
cd ..

## comparison ############################################################################
COMPdir="./table55vs45_comparison"
mkdir -p ${COMPdir}

for FILE1 in ./table45_sampleID/Summary_*.csv
do
    cp ${FILE1} ${COMPdir}/$(basename ${FILE1})
done

for FILE2 in ${OUTdir}/Summary_*.csv
do
    cp ${FILE2} ${COMPdir}/$(basename ${FILE2} .csv).checked.csv
done

##########################################################################################
# BASH DONE
