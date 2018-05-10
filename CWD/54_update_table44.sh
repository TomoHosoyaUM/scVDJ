#!/bin/bash
# usage:   bash 54_update_table44.sh
# You can edit this script to update removed data in table44, and add updated data.

echo "### $0 ###"   >> ./analysis.log.txt

##########################################################################################
OUTdir="./table54_edit44"
mkdir -p ${OUTdir}

OUTcsv=${OUTdir}/54_checked_manually.csv
echo -n > ${OUTcsv}

## edit_here_to_update ###################################################################
# comma_separated
# "Cluster,Tag"
echo "Cluster,tag" >> ${OUTcsv}                ### edit here #############################
echo "Cluster,tag" >> ${OUTcsv}                ### edit here #############################



# generate updated table #################################################################
OUTup=${OUTdir}/54_updated44_will_be_added.csv
cat ./table44_error_tag_or_IMGT/cluster_error_Tag_or_IMGT.csv | head -n 1 > ${OUTup}

N44=$(cat ./table44_error_tag_or_IMGT/cluster_error_Tag_or_IMGT.csv | grep "ROW" | wc -l )
echo "### ${N44} clusters in table44 will be removed in table55"
echo "### ${N44} clusters in table44 will be removed in table55"  >> ./analysis.log.txt

cat ${OUTcsv} | grep "ROW" > ${OUTdir}/temp.txt
N54=$(cat ${OUTdir}/temp.txt | wc -l )
echo "### ${N54} clusters updated in table54 will be added in table55"
echo "### ${N54} clusters updated in table54 will be added in table55"  >> ./analysis.log.txt

bcfile="${OUTdir}/temp.txt"
i="1"
while read -r line
do
    cluster=$(echo "${line}" | cut -d "," -f 1 )
    tag=$(echo "${line}" | cut -d "," -f 2 | sed -e 's/" "//g' )
    echo "${cluster}"
    echo "  Tag= ${tag}"
    echo "${cluster}"      >> ./analysis.log.txt
    echo "  Tag= ${tag}"   >> ./analysis.log.txt
    
    C1=$(cat ./table44_error_tag_or_IMGT/cluster_error_Tag_or_IMGT.csv | grep "${cluster}" | cut -d "," -f 1-20 )
    C22=$(cat ./table44_error_tag_or_IMGT/cluster_error_Tag_or_IMGT.csv | grep "${cluster}" | cut -d "," -f 22-26 )
    echo "${C1},${tag},${C22}" >> ${OUTup}
    
let "i++"
done < "$bcfile"

rm ${OUTdir}/temp.txt

##########################################################################################
# BASH DONE