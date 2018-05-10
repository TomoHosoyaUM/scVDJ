## usage: Rscript 26_create_consensus_table.r

system(paste0("echo @@@ 26_create_consensus_table.r >> ./analysis.log.txt"))

#################################################################################################
## date: 03/23/2018

path1="fa_consensus/"
path2="table_consensus/"

system(paste0("mkdir ",path2))

files=list.files(path1,pattern="*_ccs.fa")

for(the_file in files){
    print(the_file)
    fa=readLines(paste0(path1,the_file))
    d1=length(fa)/2
    d2=length(unlist(strsplit(fa[2],split="")))
    mat=matrix(NA,nrow=d1,ncol=d2)
    for(i in 1:d1){
        mat[i,]=unlist(strsplit(fa[i*2],split=""))
    }

    mat=mat[,apply(mat,2,function(x){any(x!="-")})] # remove all-gap positions
    fa_consensus=apply(mat,2,function(x){
        tbl=table(x) # correct the bug when there are only 2 CCSs and one of them is "-"
        tbl=tbl[c(which(names(tbl)!="-"),which(names(tbl)=="-"))]
        names(which.max(tbl)) })

    output=matrix(NA,nrow=d1,ncol=10)
    colnames(output)=c("del_pos","del_count","ins_pos","ins_count","sub_pos","sub_count","mut_pos","mut_count","length","accurate")
    ind1=fa_consensus=="-" # gap
    ind2=fa_consensus!="-" # non-gap
    for (j in 1:d1){
        t0=mat[j,]!=fa_consensus # mutation
        t1=mat[j,]=="-" & ind2 # deletion
        t2=mat[j,]!="-" & ind1 # insertion
        t3=t0 & !t1 & !t2 # substitution
        output[j,1]=paste0(which(t1),collapse=";")
        output[j,3]=paste0(which(t2),collapse=";")
        output[j,5]=paste0(which(t3),collapse=";")
        output[j,7]=paste0(which(t0),collapse=";")
        output[j,2]=sum(t1)
        output[j,4]=sum(t2)
        output[j,6]=sum(t3)
        output[j,8]=sum(t0)
    }
    output[,9]=length(fa_consensus)
    output[,10]=round(1-as.numeric(output[,"mut_count"])/length(fa_consensus),4)
    write.csv(output,file=paste0(path2,sub("fa","csv",the_file)))
}
