Brief instruction

########################################################
## Files required (tree view) ##########################
.
├── CWD
│   ├── 00_readme.txt
│   ├── 20_consensus.sh
│   ├── 21_fa_ROW_COL_DV.sh
│   ├── 23_muscle.sh
│   ├── 24_find_consensus_recursive.r
│   ├── 25_combine_identical_consensus.r
│   ├── 26_trim_rmJn.sh
│   ├── 27_create_consensus_table.r
│   ├── 40_table.sh
│   ├── 41_check_barcode_sampleID.sh
│   ├── 42_cluster.r
│   ├── 43_error_potential.sh
│   ├── 44_error_tag_or_IMGT.r
│   ├── 45_sampleID.sh
│   ├── 45a_remove_error.r
│   ├── 45b_sampleID.r
│   ├── 46_sampleID_woTg.sh
│   ├── 46a_remove_Tg.r
│   ├── 50_table.sh
│   ├── 51_manual_check.sh
│   ├── 55_sampleID_checked.sh
│   ├── 55a_remove_checked.r
│   ├── 56_sampleID_checked_woTg.sh
│   ├── 56a_remove_checked_Tg.r
│   └── IMGT.txz		# at step 3
├── barcode
│   ├── DV_seq.txt
│   ├── bc_col_rc.txt
│   ├── bc_row.txt
│   ├── bc_row_rc.txt
│   └── tag_status.csv
├── barcode_sampleID
│   ├── IDcell_sampleA.fa
│   ├── IDcell_sampleB.fa
│   └── IDcell_sampleC.fa
└── input_CCS
    ├── sampleA.fa
    ├── sampleB.fa
    └── sampleC.fa


## working directory ###################################
# cd CWD   # your preferred name

## Tools required  #####################################
# samtools     (http://www.htslib.org/)
# seqkit       (http://bioinf.shenwei.me/seqkit/)
# muscle       (http://www.drive5.com/muscle)


## 1. convert PacBio bam to fasta ######################
 samtools fasta -0 OUT.fa IN.bam  # PacBio CCS bam file

## 2. demultiplex and cluster ##########################
# usage:   bash 20_consensus.sh

## 3. IMGT #############################################
# IMGT/HighV-QUEST
# http://www.imgt.org/

## 4. create customized tables #########################
# usage:   bash 40_table.sh

## 5. (optional) manual check and edit #################
# usage:   bash 50_table.sh

########################################################



## Files required to run 40_table.sh ###################
# IMGT.txz# fa_consensus
# fa_consensus_rmJn_DJ
# fa_consensus_rmJn_VDJ
# summary_of_consensus.csv
# muscle_htmlout 




