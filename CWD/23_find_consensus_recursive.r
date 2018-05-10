## usage:   Rscript 23_find_consensus_recursive.r  20 200 3000

system(paste0("echo @@@ 23_find_consensus_recursive.r >> ./analysis.log.txt"))

#################################################################################################
# 20 is the length difference cutoff for the initial hierarchical clustering to identify small clusters
# 200 & 3000 are the cutoffs of the sequence length; <200 or >3000 are excluded

## This version defines cluster by 1.seq_length_diff 2.confidence_cutoff
## Only clusters with >= confid_cutoff and >= 2 sequences are recorded

## 1130 update:
# add top 5 minimum confidence score
# reformat the table - one cluster one row
# add random seed - otherwise the kmeans subcluster part has different results

## 180407 update:
# new, line 157, hc=hclust(dist((fa_length),method="manhattan"),method="median")
# old, line 154, hc=hclust(dist((fa_length),method="manhattan"))

set.seed(3.14)

cutoff_confidence=0.51

path1="muscle_fastaout_nowrap/"
path2="fa_consensus/"

system(paste0("mkdir ",path2))

files=list.files(path1)

args <- commandArgs(trailingOnly = TRUE)
#if(length(args)!=3){print("Usage: Rscript ind_consensus.r num_cluster num_min num_max")}

diff_length=as.numeric(args[1])

if(length(args)==3){
        num_min=as.numeric(args[2])
        num_max=as.numeric(args[3])
}else{ # if only one input, include all sequences
        num_min=0
        num_max=9999999
}

## This function cluster sequences recursively till all groups pass confid cutoff with >=2 sequences
recursive_consensus <- function(mat,list_info,cutoff,grpid,filename){
        output=list(grpid=grpid,list_info=list_info)
        tmp_consensus=consensus_seq(mat)
        if(tmp_consensus$min_confid >= cutoff){
                tmp=tmp_consensus
                grpid=grpid+1
                write_consensus(tmp,grpid,filename=filename)
                list_info$num_seq=c(list_info$num_seq,tmp$num_seq)
                list_info$seq_length=c(list_info$seq_length,tmp$seq_length)
                list_info$min_confid=c(list_info$min_confid,tmp$min_confid)
                list_info$min1=c(list_info$min1,tmp$min1)
                list_info$min2=c(list_info$min2,tmp$min2)
                list_info$min3=c(list_info$min3,tmp$min3)
                list_info$min4=c(list_info$min4,tmp$min4)
                list_info$min5=c(list_info$min5,tmp$min5)
                output=list(grpid=grpid,list_info=list_info)
        }else if(tmp_consensus$num_seq>2){ # if only 2 different sequences, drop
                tmp_grps=cluster_kmean(mat)
                for(i in 1:length(unique(tmp_grps))){
                        if(sum(tmp_grps==i)>1){
                                output=recursive_consensus(mat[tmp_grps==i,],output$list_info,cutoff,
                                        output$grpid,filename)
                        }
                }
        }
        return(output)
}

## This function excludes short/long sequences
clean_fa <- function(fa,lmin,lmax){
        ind=NULL
        if(lmin>=lmax){ print("Warning! lmax should be larger than lmin!") }
        for(i in 1:(length(fa)/2)){
                l=sum(unlist(strsplit(fa[i*2],split=""))!="-")
                if(l>=lmin & l<=lmax){ ind=c(ind,i) }
        }
        return(fa[sort(c(ind*2-1,ind*2))])
}

## This function caculate consensus information of input matrix - columns are bp positions              
consensus_seq <- function(mat){
        mat_raw=mat
        mat=mat[,apply(mat,2,function(x){any(x!="-")})] # remove all-gap positions
        #fa_consensus=apply(mat,2,function(x){ names(which.max(table(x))) }) # bug version!!
        fa_consensus=apply(mat,2,function(x){ 
            tbl=table(x) # correct the bug when there are only 2 CCSs and one of them is "-"
            tbl=tbl[c(which(names(tbl)!="-"),which(names(tbl)=="-"))] 
            names(which.max(tbl)) })
        fa_confid=apply(mat,2,function(x){ max(table(x))/length(x) })
        fa_confid=fa_confid[fa_consensus!="-"] # if a gap "-" dominates, it's an insertion and removed for now
        fa_consensus=fa_consensus[fa_consensus!="-"]
        num_seq=dim(mat)[1]
        seq_length=length(fa_consensus)
        min_confid=round(sort(fa_confid)[1],3)
        tmp=round(sort(fa_confid),3)
        min1=tmp[1];min2=tmp[2];min3=tmp[3];min4=tmp[4];min5=tmp[5]
        output=list(fa_consensus=fa_consensus,fa_confid=fa_confid,num_seq=num_seq,
                seq_length=seq_length,min_confid=min_confid,
                min1=min1,min2=min2,min3=min3,min4=min4,min5=min5,
                raw_reads=mat_raw)
        return(output)
}

## This function performs kmean cluster of input seq matrix;
## the input sequences should NOT be identical otherwise kmeans error
cluster_kmean <- function(mat,num_cluster=2){
        mat=toupper(mat)
        d1=dim(mat)[1];d2=dim(mat)[2]
        mat_bin=matrix(0,nrow=d1,ncol=d2*4)
        for(i in 1:d1){
                seq_bin=NULL
                seq_bin=c(seq_bin,mat[i,]=="A")
                seq_bin=c(seq_bin,mat[i,]=="T")
                seq_bin=c(seq_bin,mat[i,]=="G")
                seq_bin=c(seq_bin,mat[i,]=="C")
                mat_bin[i,]=as.numeric(seq_bin)
        }
        grps=kmeans(mat_bin,num_cluster)$cluster
        return(grps)
}

## This function write consensus fa and confidence csv
write_consensus <- function(list_consensus,num,filename){
        id=paste0(">",sub("\\.fa",paste0("_cluster",num),basename(filename)))
        writeLines(c(id,paste0(list_consensus$fa_consensus,collapse="")),
                con=sub("\\.fa",paste0("_cluster",num,"_consensus.fa"),filename))
        write.csv(rbind(list_consensus$fa_consensus,list_consensus$fa_confid),
                file=sub("\\.fa",paste0("_cluster",num,".csv"),filename),quote=F,row.names=F)
        # raw reads
        mat=list_consensus$raw_reads
        raw_fa=NULL
        for(i in 1:dim(mat)[1]){
                raw_fa=c(raw_fa,paste0(">",i)) # name it 1,2,3,..
                raw_fa=c(raw_fa,paste0(mat[i,],collapse=""))
        }
        writeLines(raw_fa, con=sub("\\.fa",paste0("_cluster",num,"_ccs.fa"),filename))
}

tbl_summary=NULL
for(the_file in files){
        print(the_file)
        fa=readLines(paste0(path1,the_file))
        fa=clean_fa(fa,lmin=num_min,lmax=num_max)
        if((length(fa)/2)<2) next # if only 0/1 sequence, skip
        ids=fa[seq(1,length(fa),2)]

        d1=length(fa)/2
        d2=length(unlist(strsplit(fa[2],split="")))
        mat_seq=matrix(NA,nrow=d1,ncol=d2)
        for(i in 1:d1){
                mat_seq[i,]=unlist(strsplit(fa[i*2],split=""))
        }

        # hclust by the sequence length - if two sequences are >= diff_length, they are clustered into different groups
        fa_length=unlist(lapply(as.list(fa[seq(2,length(fa),2)]),function(x){
                sum(unlist(strsplit(x,split=""))!="-")
                }))
        hc=hclust(dist((fa_length),method="manhattan"),method="median")
        grps=cutree(hc,h=diff_length)

        # process each cluster
        num_grps=length(unique(grps))
        list_info=list(num_seq=NULL,seq_length=NULL,min_confid=NULL,
                min1=NULL,min2=NULL,min3=NULL,min4=NULL,min5=NULL)
        list_final=list(grpid=0,list_info=list_info)
        for(k in 1:num_grps){
                ind=which(grps==k)
                if(length(ind)>1){ # if a cluster only has one sequence, it is excluded
                        mat=mat_seq[grps==k,]
                        # recursive consensus
                        list_final=recursive_consensus(mat,list_final$list_info,cutoff_confidence,
                                list_final$grpid,paste0(path2,the_file))
                }
        }
        info=c(sub("\\.fa","",the_file), sum(list_final$list_info$num_seq), list_final$grpid,
                paste0(list_final$list_info$num_seq,collapse=";"), paste0(list_final$list_info$seq_length,collapse=";"),
                paste0(list_final$list_info$min1,collapse=";"),
                paste0(list_final$list_info$min2,collapse=";"),
                paste0(list_final$list_info$min3,collapse=";"),
                paste0(list_final$list_info$min4,collapse=";"),
                paste0(list_final$list_info$min5,collapse=";"))
        tbl_summary=rbind(tbl_summary,info)
}

colnames(tbl_summary)=c("barcode","total_number_of_reads","number_of_clusters",
        "number_of_reads_in_each_cluster","consensus_sequence_length",
        "minimum_consensus_score1","minimum_consensus_score2","minimum_consensus_score3",
        "minimum_consensus_score4","minimum_consensus_score5")
# exclude cells without any consensus sequence
tbl_summary=tbl_summary[tbl_summary[,"total_number_of_reads"]!="0",]

#write.csv(tbl_summary,file="summary_of_consensus_sequences.csv",quote=F,row.names=F)

## reformat - one cluster one row
tbl_new=matrix(NA,nrow=0,ncol=8)
colnames(tbl_new)=c("barcode","number_of_reads","consensus_sequence_length",
        "minimum_consensus_score1","minimum_consensus_score2","minimum_consensus_score3",
        "minimum_consensus_score4","minimum_consensus_score5")
for(i in 1:dim(tbl_summary)[1]){
        l=as.numeric(tbl_summary[i,"number_of_clusters"])
        tmp=matrix(NA,nrow=l,ncol=dim(tbl_new)[2])
        for(k in 1:l){
                tmp[k,1]=paste0(tbl_summary[i,"barcode"],"_cluster",k)
        }
        tmp[,2]=unlist(strsplit(tbl_summary[i,"number_of_reads_in_each_cluster"],split=";"))
        tmp[,3]=unlist(strsplit(tbl_summary[i,"consensus_sequence_length"],split=";"))
        tmp[,4]=unlist(strsplit(tbl_summary[i,"minimum_consensus_score1"],split=";"))
        tmp[,5]=unlist(strsplit(tbl_summary[i,"minimum_consensus_score2"],split=";"))
        tmp[,6]=unlist(strsplit(tbl_summary[i,"minimum_consensus_score3"],split=";"))
        tmp[,7]=unlist(strsplit(tbl_summary[i,"minimum_consensus_score4"],split=";"))
        tmp[,8]=unlist(strsplit(tbl_summary[i,"minimum_consensus_score5"],split=";"))
        tbl_new=rbind(tbl_new,tmp)
}

write.csv(tbl_new,file="summary_of_consensus.csv",quote=F,row.names=F)














