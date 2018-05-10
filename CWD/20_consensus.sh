#!/bin/bash
# usage:   bash 20_consensus.sh

DATEstart=$(date),  UTCstart=$(date +%s)
echo "### $PWD/$0 ###"                                                     >> ./analysis.log.txt
echo "### $DATEstart ------------------------------ ### START ### $0 ###"  >> ./analysis.log.txt
echo "### $DATEstart ------------------------------ ### START ### $0 ###"  
echo " "                                                                   >> ./analysis.log.txt

######################################################################
## 1. convert bam to fasta & concatenate #############################
######################################################################
## 2. demultiplex and cluster ########################################

# need  ../input_CCS/*.fa
Nccs=$(cat ../input_CCS/*.fa | grep ">" | wc -l )
if [ ${Nccs} -gt 0 ]; then 
    echo "*** ${Nccs} CCSs found"
else
    echo "*** ERROR   ../input_CCS/*.fa does not exist, or is empty"
    echo "*** exit  $0 "
    exit
fi

bash 21_fa_ROW_COL_DV.sh   # gather sequences based on ROW-COL barcodes & DV sequences
bash 22_muscle.sh          # muscle is required; this step may take several hours
Rscript 23_find_consensus_recursive.r  20 200 3000   # find_consensus
Rscript 24_combine_identical_consensus.r             # merge identical sequences from different clusters
bash 25_trim_rmJn.sh                                 # remove adapter and extra J region

# (optional) ### calculate deletion/insertion/substitution of each raw CCS
# Rscript 26_create_consensus_table.r      


######################################################################
## 3. IMGT ###########################################################
  # IMGT/HighV-QUEST
  # http://www.imgt.org/
######################################################################
## 4. create customized tables #######################################
######################################################################
## 5. (optional) remove clusters for Tg ##############################
######################################################################
## 6. (optional) manual check and edit ###############################
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
