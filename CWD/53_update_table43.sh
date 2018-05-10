#!/bin/bash
# usage:   bash 53_update_table43.sh

echo "### $0 ###"   >> ./analysis.log.txt

##########################################################################################
OUTdir="./table53_edit43"
mkdir -p ${OUTdir}

OUTtxt=${OUTdir}/53_checked_manually.txt
echo -n > ${OUTtxt}

## edit_here_to_update ###################################################################
# "Cluster"
echo "Cluster" >> ${OUTtxt}                ### edit here #############################
echo "Cluster" >> ${OUTtxt}                ### edit here #############################



# generate updated table #################################################################
OUTup=${OUTdir}/53_updated43_will_be_removed.txt
cat ${OUTtxt} | grep "ROW" > ${OUTup}

N53=$(cat ${OUTup} | wc -l )
echo "### ${N53} clusters in table53 will be removed in table55"  >> ./analysis.log.txt
echo "### ${N53} clusters in table53 will be removed in table55"

cat ${OUTup}                                                      >> ./analysis.log.txt
cat ${OUTup}

##########################################################################################
# BASH DONE