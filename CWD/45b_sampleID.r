## usage:   Rscript 45b_sampleID.r ../../barcode_sampleID/IDcell_*.txt

system(paste0("echo @@@ 45b_sampleID.r"))
system(paste0("echo @@@ 45b_sampleID.r >> ../analysis.log.txt"))

args <- commandArgs(trailingOnly = TRUE)

#################################################################################################
## input ########################################################################################
tbl=as.matrix(read.csv("./Cluster_selected.csv",check.names=F))
tbl=tbl[which(tbl[,1]!="0"),]    ###
rownames(tbl)=tbl[,1]

tblstatus=as.matrix(read.csv("../../barcode/tag_status.csv",check.names=F))
rownames(tblstatus)=tblstatus[,"TagCombined"]

## output files #################################################################################
OUT5A=c("SampleID","GL_GL","DJ_GL","DJ_DJ","uVDJ_GL","pVDJ_GL","uVDJ_DJ","pVDJ_DJ","uVDJ_uVDJ","pVDJ_uVDJ","pVDJ_pVDJ","pVDJ_pVDJ","ID12","ID13","ID14","ID15","Partial")
OUT5B=c("SampleID","GL_GL","DJ_GL","DJ_DJ","uVDJ_GL","pVDJ_GL","uVDJ_DJ","pVDJ_DJ","uVDJ_uVDJ","pVDJ_uVDJ","pVDJ_pVDJ","ID12to15")
OUT2T=NULL
OUT2V=NULL
OUT2J=NULL

## generate tables ##############################################################################
for(filename in args){
    print(filename)
    id1=sub("IDcell_","",basename(filename))
    id2=sub(".txt","",id1)
    list_cell=readLines(filename)
    
#1. Cluster_*.csv
    ind=NULL
    for(i in list_cell){
        ind=c(ind,grep(paste0(i,"\\."),rownames(tbl)))
    }
    out1=tbl[ind,]
    write.csv(out1,file=paste0("Cluster_",id2,".csv"),row.names=F)
    
#2. Sequence_*.csv
    outVDJ<-subset(out1[,c(1:2, 16:18, 20:25)], as.numeric(out1[,"DV"]) >=3)
    outDJ<-subset(out1[,c(1:2, 20:25)], as.numeric(out1[,"DV"]) <=2)
    outVDJs<-outVDJ[order(outVDJ[,"DV"], outVDJ[,"rmJn"]),]
    outDJs<-outDJ[order(outDJ[,"Tag"], outDJ[,"rmJn"]),]
    write.csv(outVDJs,file=paste0("Sequence_VDJ_",id2,".csv"),row.names=F)
    write.csv(outDJs,file=paste0("Sequence_DJ_",id2,".csv"),row.names=F)
    
#3. CELLcluster_*.csv
    out3=matrix(NA,nrow=0,ncol=dim(tbl)[2])
    colnames(out3)=colnames(tbl)
    for(i in list_cell){
        out3=rbind(out3,NA,NA)
        out3[dim(out3)[1]-1,1]="---"
        out3[dim(out3)[1],1]=i
        ind=grep(paste0(i,"\\."),rownames(tbl))
        if(length(ind)>1){
            tmp=tbl[ind,"Vregion"]
            tmp[is.na(tmp)]="zzzzzz"
            out3=rbind(out3,tbl[ind[order(tmp)],])
        }else if(length(ind)==1){
            out3=rbind(out3,tbl[ind,])
        }
    }
    write.csv(out3,file=paste0("CELLcluster_",id2,".csv"),na="",row.names=F)
    
#4. CELL_*.csv
    z5=matrix("zzzzz", ncol=1, nrow=4) 
    head4=c("barcode","Tag1","Tag2","Tag3","Tag4","Pattern","Status","ID","TagCombined")
    out4=matrix("NA",nrow=length(list_cell),ncol=length(head4))
    rownames(out4)=list_cell
    colnames(out4)=head4
    out4[,"barcode"]=list_cell

    for(i in list_cell){
        t1=grep(paste0(i,"\\."),rownames(tbl))
        t2=as.matrix(tbl[t1,"Tag"])
        t3=rbind(t2, z5)
        t4=as.matrix(t3[order(t3[,1]), ])
        rownames(t4)=NULL
        t5=as.matrix(t4[1:4, ])
        t6=t(t5)
        t7=gsub("zzzzz", "__", t6)
        out4[i,2:5]=t7[1,1:4]

        t8=paste(t7, collapse=".")
        rownames(tblstatus)=tblstatus[,"TagCombined"]
        for(j in which(tblstatus[,"TagCombined"] %in% t8)){
            pat=tblstatus[j,5:7]
            out4[i,6:9]=tblstatus[j,5:8]
        }
    }
    write.csv(out4,file=paste0("CELL_",id2,".csv"),row.names=F)
    
#4t. count TagCombined from out4
    out4t=as.matrix(transform(tblstatus, "Num_cells"="num"))
    rownames(out4t)=out4t[, "TagCombined"]
    for(tagc in rownames(tblstatus)){
        x=subset(out4, out4[, "TagCombined"]==tagc)
        out4t[tagc,"Num_cells"]=nrow(x)
    }
    write.csv(out4t,file=paste0("Status_full_",id2,".csv"),row.names=F)

    out4s=subset(out4t, out4t[, "Num_cells"] > 0)
    write.csv(out4s,file=paste0("Status_short_",id2,".csv"),row.names=F)
    
#5. Count status from out4
    head5=c("SampleID","ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10","ID11","ID12","ID13","ID14","ID15","IDNA")
    out5a=matrix("NA",nrow=1,ncol=length(head5))
    rownames(out5a)=id2
    colnames(out5a)=head5
    out5a[1,"SampleID"]=id2
    
    st1=subset(out4, as.numeric(out4[,"ID"]) =="1")   # output is OK, but show error
    st2=subset(out4, as.numeric(out4[,"ID"]) =="2")
    st3=subset(out4, as.numeric(out4[,"ID"]) =="3")
    st4=subset(out4, as.numeric(out4[,"ID"]) =="4")
    st5=subset(out4, as.numeric(out4[,"ID"]) =="5")
    st6=subset(out4, as.numeric(out4[,"ID"]) =="6")
    st7=subset(out4, as.numeric(out4[,"ID"]) =="7")
    st8=subset(out4, as.numeric(out4[,"ID"]) =="8")
    st9=subset(out4, as.numeric(out4[,"ID"]) =="9")
    st10=subset(out4, as.numeric(out4[,"ID"]) =="10")
    st11=subset(out4, as.numeric(out4[,"ID"]) =="11")
    st12=subset(out4, as.numeric(out4[,"ID"]) =="12")
    st13=subset(out4, as.numeric(out4[,"ID"]) =="13")
    st14=subset(out4, as.numeric(out4[,"ID"]) =="14")
    st15=subset(out4, as.numeric(out4[,"ID"]) =="15")
    stNA=subset(out4, out4[,"ID"] =="NA")
    
    out5a[1,"ID1"]=nrow(st1)
    out5a[1,"ID2"]=nrow(st2)
    out5a[1,"ID3"]=nrow(st3)
    out5a[1,"ID4"]=nrow(st4)
    out5a[1,"ID5"]=nrow(st5)
    out5a[1,"ID6"]=nrow(st6)
    out5a[1,"ID7"]=nrow(st7)
    out5a[1,"ID8"]=nrow(st8)
    out5a[1,"ID9"]=nrow(st9)
    out5a[1,"ID10"]=nrow(st10)
    out5a[1,"ID11"]=nrow(st11)
    out5a[1,"ID12"]=nrow(st12)
    out5a[1,"ID13"]=nrow(st13)
    out5a[1,"ID14"]=nrow(st14)
    out5a[1,"IDNA"]=nrow(stNA)
    out5a[1,"ID15"]=nrow(st15)
    
    st10to11=sum(nrow(st10), nrow(st11))
    st12to15=sum(nrow(st12), nrow(st13), nrow(st14), nrow(st15))
    
    OUT5A=rbind(OUT5A, out5a)
    write.csv(OUT5A,file="Summary_status_full.csv",row.names=F)
    
    out5b=c(out5a[1,1:10], "ID10to11"=st10to11, "ID12to15"=st12to15)
    OUT5B=rbind(OUT5B, out5b)
    write.csv(OUT5B,file="Summary_status_short.csv",row.names=F)
    
#1t, Count Tag from out1
    head2t=c("SampleID","D1_germ_J1","D2_germ_J2","D1_x_J1","D1_x_J2","D2_x_J2","D2_x_J1","uV_x_J1","uV_x_J2","uV31_x_J1","uV31_x_J2","pV_x_J1","pV_x_J2","pV31_x_J1","pV31_x_J2")
    out2t=matrix(NA,nrow=1,ncol=length(head2t))
    colnames(out2t)=head2t
    rownames(out2t)=id2

    out2t[1,"SampleID"  ]=id2
    out2t[1,"D1_germ_J1"]=nrow(subset(out1, out1[,"Tag"] =="D1_germ_J1"))
    out2t[1,"D2_germ_J2"]=nrow(subset(out1, out1[,"Tag"] =="D2_germ_J2"))
    out2t[1,"D1_x_J1"   ]=nrow(subset(out1, out1[,"Tag"] =="D1_x_J1"))
    out2t[1,"D1_x_J2"   ]=nrow(subset(out1, out1[,"Tag"] =="D1_x_J2"))
    out2t[1,"D2_x_J2"   ]=nrow(subset(out1, out1[,"Tag"] =="D2_x_J2"))
    out2t[1,"D2_x_J1"   ]=nrow(subset(out1, out1[,"Tag"] =="D2_x_J1"))
    
    out2t[1,"uV_x_J1"   ]=nrow(subset(out1, out1[,"Tag"] =="uV_x_J1"))
    out2t[1,"uV_x_J2"   ]=nrow(subset(out1, out1[,"Tag"] =="uV_x_J2"))
    out2t[1,"uV31_x_J1" ]=nrow(subset(out1, out1[,"Tag"] =="uV31_x_J1"))
    out2t[1,"uV31_x_J2" ]=nrow(subset(out1, out1[,"Tag"] =="uV31_x_J2"))
    out2t[1,"pV_x_J1"   ]=nrow(subset(out1, out1[,"Tag"] =="pV_x_J1"))
    out2t[1,"pV_x_J2"   ]=nrow(subset(out1, out1[,"Tag"] =="pV_x_J2"))
    out2t[1,"pV31_x_J1" ]=nrow(subset(out1, out1[,"Tag"] =="pV31_x_J1"))
    out2t[1,"pV31_x_J2" ]=nrow(subset(out1, out1[,"Tag"] =="pV31_x_J2"))    
    
    OUT2T=rbind(OUT2T, out2t)
    write.csv(OUT2T,file="Summary_Tag.csv",row.names=F)
    
#1v, Count V from out1
    head2v=c("SampleID","V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12-1","V12-2","V12-3","V13-1","V13-2","V13-3","V14","V15","V16","V17","V18","V19","V20","V21","V22","V23","V24","V25","V26","V27","V28","V29","V30","V31")
    out2v=matrix(NA,nrow=1,ncol=length(head2v))
    colnames(out2v)=head2v
    rownames(out2v)=id2

    out2v[1,"SampleID"]=id2
    out2v[1,"V1"   ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV1"))
    out2v[1,"V2"   ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV2"))
    out2v[1,"V3"   ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV3"))
    out2v[1,"V4"   ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV4"))
    out2v[1,"V5"   ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV5"))
    out2v[1,"V6"   ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV6"))
    out2v[1,"V7"   ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV7"))
    out2v[1,"V8"   ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV8"))
    out2v[1,"V9"   ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV9"))
    out2v[1,"V10"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV10"))
    
    out2v[1,"V11"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV11"))
    out2v[1,"V12-1"]=nrow(subset(out1, out1[,"Vregion"] =="TRBV12-1"))
    out2v[1,"V12-2"]=nrow(subset(out1, out1[,"Vregion"] =="TRBV12-2"))
    out2v[1,"V12-3"]=nrow(subset(out1, out1[,"Vregion"] =="TRBV12-3"))
    out2v[1,"V13-1"]=nrow(subset(out1, out1[,"Vregion"] =="TRBV13-1"))
    out2v[1,"V13-2"]=nrow(subset(out1, out1[,"Vregion"] =="TRBV13-2"))
    out2v[1,"V13-3"]=nrow(subset(out1, out1[,"Vregion"] =="TRBV13-3"))
    out2v[1,"V14"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV14"))
    out2v[1,"V15"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV15"))
    out2v[1,"V16"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV16"))
    out2v[1,"V17"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV17"))
    out2v[1,"V18"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV18"))
    out2v[1,"V19"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV19"))
    out2v[1,"V20"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV20"))
    
    out2v[1,"V21"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV21"))
    out2v[1,"V22"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV22"))
    out2v[1,"V23"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV23"))
    out2v[1,"V24"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV24"))
    out2v[1,"V25"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV25"))
    out2v[1,"V26"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV26"))
    out2v[1,"V27"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV27"))
    out2v[1,"V28"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV28"))
    out2v[1,"V29"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV29"))
    out2v[1,"V30"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV30"))
    out2v[1,"V31"  ]=nrow(subset(out1, out1[,"Vregion"] =="TRBV31"))

    OUT2V=rbind(OUT2V, out2v)
    write.csv(OUT2V,file="Summary_region_V_in_VDJ.csv",row.names=F)

#1v, Count J from out1
    head2j=c("SampleID","J1-1","J1-2","J1-3","J1-4","J1-5","J1-6","J1-7","J2-1","J2-2","J2-3","J2-4","J2-5","J2-6","J2-7")
    out2j=matrix(NA,nrow=1,ncol=length(head2j))
    colnames(out2j)=head2j
    rownames(out2j)=id2
    
    out2j[1,"SampleID"]=id2
    out2j[1,"J1-1"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ1-1"))
    out2j[1,"J1-2"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ1-2"))
    out2j[1,"J1-3"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ1-3"))
    out2j[1,"J1-4"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ1-4"))
    out2j[1,"J1-5"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ1-5"))
    out2j[1,"J1-6"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ1-6"))
    out2j[1,"J1-7"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ1-7"))

    out2j[1,"J2-1"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ2-1"))
    out2j[1,"J2-2"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ2-2"))
    out2j[1,"J2-3"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ2-3"))
    out2j[1,"J2-4"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ2-4"))
    out2j[1,"J2-5"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ2-5"))
    out2j[1,"J2-6"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ2-6"))
    out2j[1,"J2-7"]=nrow(subset(out1, out1[,"Jregion"] =="TRBJ2-7"))

    OUT2J=rbind(OUT2J, out2j)
    write.csv(OUT2J,file="Summary_region_J_in_VDJ.csv",row.names=F)  
    
}
