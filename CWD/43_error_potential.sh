#!/bin/bash
# usage:   bash 43_error_potential.sh

echo "### $0 ###"   >> ./analysis.log.txt
echo "### $0 ###"  

#################################################################################################
OUTdir=./table43_error_potential
mkdir -p ${OUTdir}

OUTmulti=${OUTdir}/multi_cluster_ROW.COL.DV.txt
    echo -n > $OUTmulti

OUTcsv=${OUTdir}/multi_cluster.csv
cat ./table42_Cluster.csv | head -n 1 > ${OUTcsv}

CHECK=${OUTdir}/need_manual_check_cluster.csv
cat ./table42_Cluster.csv | head -n 1 > ${CHECK}

echo -n > ${OUTdir}/need_manual_check_cluster.txt


# only when multiple cluster exist ##############################################################
C2=$(find ./fa_consensus/*cluster2_consensus.fa | wc -l ) 
if [ ${C2} -gt 0 ]; then
    echo "### ./fa_consensus/*cluster2_consensus.fa  exist"
else
    echo "### No *cluster2_consensus.fa" >> ./analysis.log.txt
    echo "### No *cluster2_consensus.fa"
    echo "### skip $0 "
    exit
fi

for FILE in ./fa_consensus/*cluster2_consensus.fa
do
    RCDV=$(basename $FILE _cluster2_consensus.fa)
# generate list for barcode with 2 clusters
    echo "${RCDV}" >> $OUTmulti
# copy Cluster info
    echo "---" >> ${OUTcsv}
    cat ./table42_Cluster.csv | grep "${RCDV}_" >> ${OUTcsv}

# using same DV and J
    vj=$(cat ./table42_Cluster.csv | grep "${RCDV}_" | tail -n 1 | cut -d ',' -f 22,23 )
    dup=$(cat ./table42_Cluster.csv | grep "${RCDV}_" | grep "${vj}" | wc -l )
    if [${dup} -eq 2 ]; then 
        echo "---" >> ${CHECK}
        cat ./table42_Cluster.csv | grep "${RCDV}_" >> ${CHECK}
    fi
# >3 clusters
    N=$(cat ./table42_Cluster.csv | grep "${RCDV}_" | wc -l )
    if [${N} -gt 2 ]; then 
        echo "---" >> ${CHECK}
        cat ./table42_Cluster.csv | grep "${RCDV}_" >> ${CHECK}
    fi
done

# need_manual_check
cat ${CHECK} | cut -d ',' -f 2 | sed -e 1d | grep -v "^---" | sed -e 's/"//g' > ${OUTdir}/need_manual_check_cluster.txt
cat ${OUTdir}/need_manual_check_cluster.txt | sed -e 's/_cluster.*//g' | sort | uniq > ${OUTdir}/need_manual_check_barcode.txt

bcfile="${OUTdir}/need_manual_check_barcode.txt"
i="1"
while read -r line
do
    echo "$line"
# copy consensus.fa
    cat ./fa_consensus/${line}_*_consensus.fa >  ${OUTdir}/${line}_consensus.fa
# muscle consensus
    muscle -in ${OUTdir}/${line}_consensus.fa -htmlout ${OUTdir}/${line}_consensus.html
    
let "i++"
done < "$bcfile"

# log
Ntotal=$(cat ${OUTcsv} | grep "_cluster" | wc -l )
Ncheck=$(cat ${CHECK}  | grep "_cluster" | wc -l )

echo "### multiple clusters were generated in ${C2} wells"    >> ./analysis.log.txt
echo "### ${Ntotal} clusters in total "                       >> ./analysis.log.txt
echo "### ${Ncheck} clusters need to be checked manually "    >> ./analysis.log.txt

echo "### multiple clusters were generated in ${C2} wells"
echo "### ${Ntotal} clusters in total"
echo "### ${Ncheck} clusters need to be checked manually, removing for now"

#################################################################################################
# BASH DONE
