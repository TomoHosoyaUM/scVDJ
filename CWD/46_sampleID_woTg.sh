#!/bin/bash
# usage:   bash 46_sampleID_woTg.sh

echo "### $0 ###"   >> ./analysis.log.txt
echo "### $0 ###" 

##########################################################################################
OUTdir="./table46_sampleID_woTg"
mkdir -p ${OUTdir}

Rscript 46a_remove_Tg.r

cd ${OUTdir}
Rscript ../45b_sampleID.r ../../barcode_sampleID/IDcell_*.txt
cd ..

## comparison ############################################################################
COMPdir="./table46vs45_comparison"
mkdir -p ${COMPdir}

for FILE1 in ./table45_sampleID/Summary_*.csv
do
    cp ${FILE1} ${COMPdir}/$(basename ${FILE1})
done

for FILE2 in ${OUTdir}/Summary_*.csv
do
    cp ${FILE2} ${COMPdir}/$(basename ${FILE2} .csv).woTg.csv
done

##########################################################################################
# BASH DONE
