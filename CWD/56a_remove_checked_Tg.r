## usage:   Rscript 56a_remove_checked_Tg.r

system(paste0("echo @@@ 56a_remove_checked_Tg.r >> ./analysis.log.txt"))
system(paste0("echo @@@ 56a_remove_checked_Tg.r"))

#################################################################################################
## output directory #############################################################################
path1="./table56_sampleID_checked_woTg/"
system(paste0("mkdir -p ",path1))

## for all clusters #############################################################################
tbl=as.matrix(read.csv("./table42_Cluster.csv",check.names=F))
tbl=tbl[which(tbl[,2]!="0"),]    ###
rownames(tbl)=tbl[,2]
tbl=tbl[,-1]

Nall=nrow(tbl)
print(Nall)
system(paste0("echo ",Nall," total clusters in ./table42_Cluster.csv >> ./analysis.log.txt"))

## remove table53_edit43 ########################################################################
list_remove=readLines("./table53_edit43/53_updated43_will_be_removed.txt")
tblerror53=tbl[(rownames(tbl) %in% list_remove),]
write.csv(tblerror53,file=paste0(path1,"53_checked43_Cluster_removed.csv"),row.names=F)

tbl=tbl[!(rownames(tbl) %in% list_remove),]

Ner53=nrow(tblerror53)
print(Ner53)
system(paste0("echo ",Ner53," clusters in table53_updated43 were removed >> ./analysis.log.txt"))

## remove tag error in table44 ##################################################################
tblerror44=as.matrix(read.csv("./table44_error_tag_or_IMGT/Cluster_error_Tag_or_IMGT.csv",check.names=F))
tblerror44=tbl[(rownames(tbl) %in% tblerror44[,2]),]
write.csv(tblerror44,file=paste0(path1,"44_error_Cluster_removed.csv"),row.names=F)

tbl=tbl[!(rownames(tbl) %in% rownames(tblerror44)),]

Ner44=nrow(tblerror44)
print(Ner44)
system(paste0("echo ",Ner44," clusters in table44 were removed >> ./analysis.log.txt"))

## add table54_edit44 ###########################################################################
tbladd54=as.matrix(read.csv("./table54_edit44/54_updated44_will_be_added.csv",check.names=F))
rownames(tbladd54)=tbladd54[,2]
tbladd54=tbladd54[,-1]
write.csv(tbladd54,file=paste0(path1,"54_updated44_Cluster_added.csv"),row.names=F)

tbl=rbind(tbladd54, tbl)
tbl=tbl[which(tbl[,1]!="0"),]    ###
rownames(tbl)=tbl[,1]

Nadd54=nrow(tbladd54)
print(Nadd54)
system(paste0("echo ",Nadd54," clusters in table54_edit44 were updated >> ./analysis.log.txt"))

## remove Tg ###################################################################################
# part of nucleotide sequence for Vb8 transgene from TRBV13-2 primer to J1-1 (221 bp)
Tg="gcatgggctgaggctgatccattattcatatggtgctggcagcactgagaaaggagatatccctgatggatacaaggcctccagaccaagccaagagaacttctccctcattctggagttggctaccccctctcagacatcagtgtacttctgtgccagcggctccgggacaacaaacacagaagtcttctttggtaaaggaaccagactcacagttgtag"

tblTg=subset(tbl, tbl[, "nt"] ==Tg)
tbl=tbl[!(rownames(tbl) %in% tblTg[, 1]),]

Ntg=nrow(tblTg)
print(Ntg)
system(paste0("echo ",Ntg," Tg clusters removed >> ./analysis.log.txt"))

## final tbl write.csv ##########################################################################
write.csv(tbl,file=paste0(path1,"56_Cluster_selected.csv"),row.names=F)
write.csv(tbl,file=paste0(path1,"Cluster_selected.csv"),row.names=F)

N=nrow(tbl)
print(N)
system(paste0("echo ",N," clusters selected in table56 >> ./analysis.log.txt"))

#################################################################################################




