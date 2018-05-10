#!/bin/bash
# usage:   bash 50_table.sh

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
######################################################################
## 5. (optional) manual check and edit ###############################

bash 53_update_table43.sh           # edit this script
bash 54_update_table44.sh           # edit this script
bash 55_sampleID_checked.sh         # table for each sample

## (optional) remove Tg ## 
# bash 56_sampleID_checked_woTg.sh    # table for each sample

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
