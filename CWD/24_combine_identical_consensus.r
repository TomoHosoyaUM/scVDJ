## usage:   Rscript 24_combine_identical_consensus.r

system(paste0("echo @@@ 24_combine_identical_consensus.r >> ./analysis.log.txt"))

#################################################################################################
## date: 03/24/2018
## This is a patch for find_consensus_recursive.r, to combine identical consensus sequences.
## find_consensus_recursive.r has two-step cluster - the 1st hcluster based on length differences is necessary
## Otherwise small clusters can not be identified at the 2nd kmeans cluster - we need to pre-assign a cluter number to kmeans (default=2)
## This patch just check multiple cluster cases. If any identical consensus, combine them, re-write corresponding files and modify the summary csv.


path1="fa_consensus/"

tbl=read.csv("summary_of_consensus.csv",stringsAsFactors=F)
tbl_new=tbl
ids=sub("_cluster.*","",tbl[,"barcode"])

ids_multi=unique(ids[duplicated(ids)])

ind_del=NULL # index of rows to be deleted
for(i in 1:length(ids_multi)){
    ind=grep(paste0(ids_multi[i],"_cluster"),tbl[,"barcode"]) # multiple clusters
    if(any(duplicated(tbl[ind,"consensus_sequence_length"]))){ # if any same length clusters, further read fa and check
        fa=NULL
        for(j in ind){
            fa=c(fa,readLines(paste0(path1,tbl[j,"barcode"],"_consensus.fa"))[2])
        }
        tbl_seq=table(fa) # check if any duplicated sequences
        if(any(tbl_seq!=1)){
            ind_new=which(tbl_seq!=1) # duplicated consensus
            for(j in ind_new){
                ind_grps=(fa==names(tbl_seq)[j]) # grps/clusters have the same consensus
                ind_del=c(ind_del,ind[ind_grps][-1]) # to be deleted
                ccs=NULL
                for(k in ind[ind_grps]){
                    ccs=c(ccs,readLines(paste0(path1,tbl[k,"barcode"],"_ccs.fa")))
                }
                d1=length(ccs)/2
                d2=nchar(ccs[2])
                ccs[c(1:d1)*2-1]=paste0(">",1:d1) # rename them to >1, >2, ..
                mat=matrix(NA,nrow=d1,ncol=d2)
                for(k in 1:d1){
                    mat[k,]=unlist(strsplit(ccs[k*2],split=""))
                }
                mat=mat[,apply(mat,2,function(x){any(x!="-")})] # remove all gap positions
                fa_consensus=apply(mat,2,function(x){
                    tmp=table(x) # correct the bug when there are only 2 CCSs and one of them is "-"
                    tmp=tmp[c(which(names(tmp)!="-"),which(names(tmp)=="-"))]
                    names(which.max(tmp)) })
                fa_confid=apply(mat,2,function(x){ max(table(x))/length(x) })
                fa_confid=fa_confid[fa_consensus!="-"] # if a gap "-" dominates, it's an insertion and removed for now
                fa_consensus=fa_consensus[fa_consensus!="-"]
                min_confid=round(sort(fa_confid),3)
                tbl[ind[ind_grps][1],"number_of_reads"]=d1
                tbl[ind[ind_grps][1],"consensus_sequence_length"]=length(fa_consensus)
                tbl[ind[ind_grps][1],paste0("minimum_consensus_score",1:5)]=min_confid[1:5]
                # write new ccs/consensue
                tmp_id=tbl[ind[ind_grps][1],"barcode"]
                print(tmp_id)
                writeLines(c(paste0(">",tmp_id),paste0(fa_consensus,collapse="")),
                    con=paste0(path1,tmp_id,"_consensus.fa"))
                write.csv(rbind(fa_consensus,fa_confid),
                    file=paste0(path1,tmp_id,".csv"),quote=F,row.names=F)
                writeLines(ccs,con=paste0(path1,tmp_id,"_ccs.fa"))
            }
        }
    }
}

for(i in ind_del){ # remove redundant files & delete lines
    tmp_id=tbl[i,"barcode"]
    system(paste0("rm ",path1,tmp_id,"*"))
}

if(!is.null(ind_del)){
    tbl=tbl[-ind_del,]
}

write.csv(tbl,file="summary_of_consensus.csv",quote=F,row.names=F)

