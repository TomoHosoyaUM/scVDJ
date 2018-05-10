#!/bin/bash
# usage:  bash 25_trim_rmJn.sh

DATEstart=$(date),  UTCstart=$(date +%s)
echo " "                                                                   >> ./analysis.log.txt
echo "### $0 ###"                                                          >> ./analysis.log.txt
echo "### $DATEstart ------------------------------ ### START ### $0 ###"  >> ./analysis.log.txt
echo "### $DATEstart ------------------------------ ### START ### $0 ###"  

#################################################################################################
# remove barcode-adapter and extra J region sequence

# remove barcode-AF and barcode-AR3
mkdir -p fa_consensus_trim

for FILE in ./fa_consensus/*_consensus.fa
do
    OUT=./fa_consensus_trim/$(basename ${FILE} .fa).trim.fa

    # seq_ID
    echo ">$(basename ${FILE} .fa).trim" > ${OUT}
    # remove barcode-AF (32bp) and barcode-AR3       (CCTGTGTGAAATTGTTATCCGC = AR3 reverse complement)
    cat ${FILE} | grep -v '>' |  cut -c 33- | sed -e "s/CCTGTGTGAAATTGTTATCCGC.*//g" >> $OUT
done

# remove extra Jn genomic DNA
# V primers (DV3 through DV34)

mkdir -p ./fa_consensus_rmJn_VDJ

for i in `seq 1 34`
do
    echo DV${i}
    
    for FILE in ./fa_consensus_trim/ROW*.DV${i}_cluster*.fa
    do
        echo -n > temp.fa
        echo ">$(basename ${FILE} .trim.fa).J1-7rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/TAAGACAGAATCCTTAGGTATAAGGTAAGA.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J1-6rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTATGGGGGCTCCATTTCTGACTGGAGGGG.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J1-5rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTAAACTATGGGACCAAACTGGTGGGACCA.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J1-4rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTATGTAAAAGATTTCTTTCTCGGGAGGGT.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J1-3rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTAAGTTAGGGCCAAATGGCTGGGTACTGG.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J1-2rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTAAGGCCTGAGGGTCTTTGGGTGTGGGAT.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J1-1rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTAAGATATCTTTCAGGTAAATTTCCAGGT.*//g" >> temp.fa

        echo ">$(basename ${FILE} .trim.fa).J2-7rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTAAGATTCACATCTCTCGCTTCCACCCAA.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J2-6rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/CTTCTTGGCAACTGCAGCGGGGAGTTCTGG.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J2-5rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTGAGCTGGGGCCCCACGTGCGCGTTCTCA.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J2-4rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTAAGCTGGGGTATAGTTTTTGTGTTGGGT.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J2-3rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTAAGTTGGGAGCTAGTAATGAAGGGGAGG.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J2-2rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTAAGCAGGCAGCTGGGGGTACCATGGGAG.*//g" >> temp.fa
        echo ">$(basename ${FILE} .trim.fa).J2-1rm" >> temp.fa
        cat ${FILE} | grep -v '>' | sed -e "s/GTAAGAAGGCAGAGGCCATACAGGTGGGAG.*//g" >> temp.fa

        cat temp.fa | seqkit sort -w 0 -l -2 | head -n 2 > ./fa_consensus_rmJn_VDJ/$(basename ${FILE} .trim.fa).rmJn.fa

    done
done

rm temp.fa

# move DV1 and DV2 from fa_consensus_rmJn_VDJ
mkdir -p ./fa_consensus_rmJn_DJ
for FILE1 in ./fa_consensus_rmJn_VDJ/*DV1_*_consensus.rmJn.fa
do
    mv ${FILE1} ./fa_consensus_rmJn_DJ/$(basename ${FILE1}) 
done

for FILE2 in ./fa_consensus_rmJn_VDJ/*DV2_*_consensus.rmJn.fa
do
    mv ${FILE2} ./fa_consensus_rmJn_DJ/$(basename ${FILE2}) 
done

# move empty from fa_consensus_rmJn_VDJ
mkdir -p ./fa_consensus_rmJn__empty
for EMPTY in ./fa_consensus_rmJn_VDJ/*.fa.rmJn.fa
do 
	mv ${EMPTY} ./fa_consensus_rmJn__empty/$(basename ${EMPTY}) 
done

# combine .rmJn.fa for IMGT 
cat ./fa_consensus_rmJn_VDJ/*.fa > IMGT.fa

# check number of sequence 
# if Lvdj and Limgt are different, possible error in the FASTA file for IMGT, re-run this script
ccs=$(cat ../input_CCS/*.fa | grep ">" | wc -l )
barcoded=$(cat ./fa_ROW_COL/*.fa | grep ">" | wc -l )
freq=$(echo " 100 * ${barcoded} / ${ccs} " | bc )
cluster=$(find ./fa_consensus/*_consensus.fa | wc -l )
djclus=$(find ./fa_consensus_rmJn_DJ/*.rmJn.fa | wc -l )
vdjclus=$(find ./fa_consensus_rmJn_VDJ/*.rmJn.fa | wc -l )
imgtfa=$(cat ./IMGT.fa | grep ">" | wc -l)

echo "### ${ccs}  CCSs found in ../input_CCS/*.fa"                                                 >> ./analysis.log.txt
echo "### ${barcoded}  CCSs have barcode-AF and barcode-AR3  flanking each end "        >> ./analysis.log.txt
echo "###             ${freq}% of total CCSs"                                           >> ./analysis.log.txt
echo "### ${cluster}  total clusters generated in ./fa_consensus"                       >> ./analysis.log.txt
echo "### ${djclus}  DJ clusters are  in ./fa_consensus_rmJn_DJ"                        >> ./analysis.log.txt
echo "### ${vdjclus}  VDJ clusters are in ./fa_consensus_rmJn_VDJ"                      >> ./analysis.log.txt
echo "### ${imgtfa}  VDJ sequences (rmJn) are in ./IMGT.fa"                             >> ./analysis.log.txt

echo "### ${ccs}  CCSs found in ../input_CCS/*.fa"     
echo "### ${barcoded}  CCSs have barcode-AF and barcode-AR3  flanking each end " 
echo "###             ${freq}% of total CCSs"     
echo "### ${cluster}  total clusters generated in ./fa_consensus"  
echo "### ${djclus}  DJ clusters are  in ./fa_consensus_rmJn_DJ"    
echo "### ${vdjclus}  VDJ clusters are in ./fa_consensus_rmJn_VDJ"   
echo "### ${imgtfa}  VDJ sequences (rmJn) are in ./IMGT.fa" 

if [ ${vdjclus} -eq ${imgtfa} ]; then
    echo "### ./IMGT.fa is ready for IMGT/HighV-QUEST analysis" >> ./analysis.log.txt
    echo "### ./IMGT.fa is ready for IMGT/HighV-QUEST analysis"
else
    echo "### Maybe ERROR in IMGT.fa, re-run [ 27_fa_consensus_trim_rmJn.sh ]" >> ./analysis.log.txt
    echo "### Maybe ERROR in IMGT.fa, re-run [ 27_fa_consensus_trim_rmJn.sh ]"
fi

echo "total_CCSs,CCSs_ROW_COL,clusters,D-to-J,V-to-DJ" > ./table25_stats.csv
echo "${ccs},${barcoded},${cluster},${djclus},${vdjclus}" >> ./table25_stats.csv

#################################################################################################
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


# Trb       J region sequence                                       immediate 3' sequence

# Trbj1-7 = CCTGTGTTGGATGACCATGGTCTTGGAAAGGAACTTAGGTATAAGA          TAAGACAGAATCCTTAGGTATAAGGTAAGA
# Trbj1-6 = TTCCTATAATTCGCCCCTCTACTTTGCGGCAGGCACCCGGCTCACTGTGACAG   GTATGGGGGCTCCATTTCTGACTGGAGGGG
# Trbj1-5 = TAACAACCAGGCTCCGCTTTTTGGAGAGGGGACTCGACTCTCTGTTCTAG      GTAAACTATGGGACCAAACTGGTGGGACCA
# Trbj1-4 = TTTCCAACGAAAGATTATTTTTCGGTCATGGAACCAAGCTGTCTGTCCTGG     GTATGTAAAAGATTTCTTTCTCGGGAGGGT
# Trbj1-3 = TTCTGGAAATACGCTCTATTTTGGAGAAGGAAGCCGGCTCATTGTTGTAG      GTAAGTTAGGGCCAAATGGCTGGGTACTGG
# Trbj1-2 = CAAACTCCGACTACACCTTCGGCTCAGGGACCAGGCTTTTGGTAATAG        GTAAGGCCTGAGGGTCTTTGGGTGTGGGAT
# Trbj1-1 = CAAACACAGAAGTCTTCTTTGGTAAAGGAACCAGACTCACAGTTGTAG        GTAAGATATCTTTCAGGTAAATTTCCAGGT

# Trbj2-7 = CTCCTATGAACAGTACTTCGGTCCCGGCACCAGGCTCACGGTTTTAG         GTAAGATTCACATCTCTCGCTTCCACCCAA
# Trbj2-6 = CAGCCCTTGCCCTGACTGATTGGCAGCCGATTGAACAGCCTATGCGAG        CTTCTTGGCAACTGCAGCGGGGAGTTCTGG
# Trbj2-5 = AACCAAGACACCCAGTACTTTGGGCCAGGCACTCGGCTCCTCGTGTTAG       GTGAGCTGGGGCCCCACGTGCGCGTTCTCA
# Trbj2-4 = AGTCAAAACACCTTGTACTTTGGTGCGGGCACCCGACTATCGGTGCTAG       GTAAGCTGGGGTATAGTTTTTGTGTTGGGT
# Trbj2-3 = AGTGCAGAAACGCTGTATTTTGGCTCAGGAACCAGACTGACTGTTCTCG       GTAAGTTGGGAGCTAGTAATGAAGGGGAGG
# Trbj2-2 = CAAACACCGGGCAGCTCTACTTTGGTGAAGGCTCAAAGCTGACAGTGCTGG     GTAAGCAGGCAGCTGGGGGTACCATGGGAG
# Trbj2-1 = TAACTATGCTGAGCAGTTCTTCGGACCAGGGACACGACTCACCGTCCTAG      GTAAGAAGGCAGAGGCCATACAGGTGGGAG
