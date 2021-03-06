## usage:   Rscript 46a_remove_Tg.r

system(paste0("echo @@@ 46a_remove_Tg.r >> ./analysis.log.txt"))
system(paste0("echo @@@ 46a_remove_Tg.r"))

args <- commandArgs(trailingOnly = TRUE)

#################################################################################################
## output directory #############################################################################
path1="./table46_sampleID_woTg/"
system(paste0("mkdir -p ",path1))

## table for clusters ###########################################################################
tbl=as.matrix(read.csv("./table42_Cluster.csv",check.names=F))
tbl=tbl[which(tbl[,2]!="0"),]    ###
rownames(tbl)=tbl[,2]
tbl=tbl[,-1]
Nall=nrow(tbl)
print(Nall)
system(paste0("echo ",Nall," total clusters in ./table42_Cluster.csv >> ./analysis.log.txt"))

## remove error in table43_ #####################################################################
tblerror43=as.matrix(read.csv("./table43_error_potential/need_manual_check_cluster.csv",check.names=F))
tblerror43=tbl[(rownames(tbl) %in% tblerror43[,2]),]
write.csv(tblerror43,file=paste0(path1,"43_error_Cluster_removed.csv"),row.names=F)

tbl=tbl[!(rownames(tbl) %in% rownames(tblerror43)),]

Ner43=nrow(tblerror43)
print(Ner43)
system(paste0("echo ",Ner43," clusters in table43 removed >> ./analysis.log.txt"))

## remove error in table44_ #####################################################################
tblerror44=as.matrix(read.csv("./table44_error_tag_or_IMGT/Cluster_error_Tag_or_IMGT.csv",check.names=F))
tblerror44=tbl[(rownames(tbl) %in% tblerror44[,2]),]
write.csv(tblerror44,file=paste0(path1,"44_error_Cluster_removed.csv"),row.names=F)

tbl=tbl[!(rownames(tbl) %in% rownames(tblerror44)),]

Ner44=nrow(tblerror44)
print(Ner44)
system(paste0("echo ",Ner44," clusters in table44 removed >> ./analysis.log.txt"))

## remove Tg ###################################################################################
# part of nucleotide sequence for Vb8 transgene from TRBV13-2 primer to J1-1 (221 bp)
Tg="gcatgggctgaggctgatccattattcatatggtgctggcagcactgagaaaggagatatccctgatggatacaaggcctccagaccaagccaagagaacttctccctcattctggagttggctaccccctctcagacatcagtgtacttctgtgccagcggctccgggacaacaaacacagaagtcttctttggtaaaggaaccagactcacagttgtag"

tblTg=subset(tbl, tbl[, "nt"] ==Tg)
tbl=tbl[!(rownames(tbl) %in% tblTg[, 1]),]

Ntg=nrow(tblTg)
print(Ntg)
system(paste0("echo ",Ntg," Tg clusters removed >> ./analysis.log.txt"))

## final tbl write.csv ##########################################################################
write.csv(tbl,file=paste0(path1,"46_Cluster_selected.csv"),row.names=F)
write.csv(tbl,file=paste0(path1,"Cluster_selected.csv"),row.names=F)

N=nrow(tbl)
print(N)
system(paste0("echo ",N," clusters selected in table46 >> ./analysis.log.txt"))

#################################################################################################







