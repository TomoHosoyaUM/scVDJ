#!/bin/bash
# usage:   bash 56_sampleID_checked_woTg.sh

echo "### $0 ###"   >> ./analysis.log.txt
echo "### $0 ###" 

##########################################################################################
OUTdir="./table56_sampleID_checked_woTg"
mkdir -p ${OUTdir}

Rscript 56a_remove_checked_Tg.r

cd ${OUTdir}
Rscript ../45b_sampleID.r ../../barcode_sampleID/IDcell_*.txt
cd ..

## comparison ############################################################################
COMPdir="./table56vs46_comparison"
mkdir -p ${COMPdir}

for FILE1 in ./table46_sampleID_woTg/Summary_*.csv
do
    cp ${FILE1} ${COMPdir}/$(basename ${FILE1} .csv).woTg.csv
done

for FILE2 in ${OUTdir}/Summary_*.csv
do
    cp ${FILE2} ${COMPdir}/$(basename ${FILE2} .csv).woTg.checked.csv
done

##########################################################################################
# BASH DONE
