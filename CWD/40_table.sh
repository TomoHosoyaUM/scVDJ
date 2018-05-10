#!/bin/bash
# usage:   bash 40_table.sh

DATEstart=$(date),  UTCstart=$(date +%s)
echo " "                                                  >> ./analysis.log.txt
echo "### $PWD/$0 ###"                                                     >> ./analysis.log.txt
echo "### $DATEstart ------------------------------ ### START ### $0 ###"  >> ./analysis.log.txt
echo "### $DATEstart ------------------------------ ### START ### $0 ###"  

######################################################################
## 1. convert bam to fasta & concatenate #############################
######################################################################
## 2. demultiplex and cluster ########################################
######################################################################
## 3. IMGT ###########################################################
######################################################################
## 4. create customized tables #######################################

# need       ./*.txz (IMGT outout)  in the working directory (same with IMGT.fa)
COUNT=$(find ./*.txz | wc -l ) 
if [ $COUNT -eq 1 ]; then
    mkdir -p ./IMGT
    cd IMGT
    tar -xf ../*.txz
    cd ..
else
    echo "0 or 2 *.txz file exist (IMGT outout), place only 1"
    echo "### exit  $0 "
    exit
fi

# need  ../barcode_sampleID/IDcell_*.txt   (list of barcode ROW.COL.)
N=$(cat ../barcode_sampleID/IDcell_*.txt | wc -l )
if [ $N -gt 0  ]; then 
    bash 41_check_barcode_sampleID.sh
else
    echo "### ERROR   No ../barcode_sampleID/cell_*.txt  file is found"
    echo "### exit  $0 "
    exit
fi

Rscript 42_cluster.r              # table for all cluster info
bash 43_error_potential.sh        # list error
Rscript 44_error_tag_or_IMGT.r    # list error
bash 45_sampleID.sh               # table for each sample

## (optional) remove Tg ## 
bash 46_sampleID_woTg.sh          # table for each sample


######################################################################
## 5. (optional) manual check and edit ###############################
######################################################################
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
