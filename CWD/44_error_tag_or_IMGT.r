### usage:   Rscript 44_error_tag_or_IMGT.r

system(paste0("echo @@@ 44_error_tag_or_IMGT.r >> ./analysis.log.txt"))
system(paste0("echo @@@ 44_error_tag_or_IMGT.r"))

#################################################################################################
## output directory #############################################################################
path1="./table44_error_tag_or_IMGT/"
system(paste0("mkdir -p ",path1))

## input files ##################################################################################
path2="fa_consensus_rmJn_VDJ/"
    filesvdj=list.files(path2)
    clusterVDJ=sub("_consensus.rmJn.fa", "", filesvdj)
tblimgt=as.matrix(read.delim("./IMGT/4_IMGT-gapped-AA-sequences.txt",check.names=F))
    idIMGTs=sub("\\_consensus.*", "", tblimgt[,2])
tblCluster=read.csv("./table42_Cluster.csv")
    rownames(tblCluster)=tblCluster[,2]
    
## duplicate in IMGT output #####################################################################
idDup=idIMGTs[duplicated(idIMGTs)]
tblErrorDup=tblCluster[idDup, ]
if (nrow(tblErrorDup) > 0){
    tblErrorDup[,1]="duplicate_in_IMGT_output"
}

## missing in IMGT output #######################################################################
idMiss<-setdiff(clusterVDJ, idIMGTs)
tblErrorMiss=tblCluster[idMiss, ]
if (nrow(tblErrorMiss) > 0){
    tblErrorMiss[,1]="missing_in_IMGT_output"
}

## error convert to Tag #########################################################################
tblErrorTag=subset(tblCluster, tblCluster[, "Tag"] =="__")
if (nrow(tblErrorTag) > 0){
    tblErrorTag[,1]="No_Tag"
}

## prepare for manual check ######################################################################
tblErrorAll=rbind(tblErrorDup, tblErrorMiss, tblErrorTag)
write.csv(tblErrorAll,file=paste0(path1,"cluster_error_Tag_or_IMGT.csv"),row.names=F)

write(rownames(tblErrorAll), file=paste0(path1,"cluster_error_Tag_or_IMGT.txt"))

## copy fasta files ##############################################################################
for(id in rownames(tblErrorAll)){
    print(id)
    fastafile=(paste0("./fa_consensus/",id,"_consensus.fa"))

    system(paste0("cp ",fastafile," ",path1,id,"_consensus.fa"))
}



