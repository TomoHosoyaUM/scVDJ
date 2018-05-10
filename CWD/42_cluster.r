## usage:  Rscript 42_cluster.r

system(paste0("echo @@@ 42_table.r >> ./analysis.log.txt"))
system(paste0("echo @@@ 42_table.r"))

#################################################################################################
## input files ##################################################################################
ccs=as.matrix(read.csv("summary_of_consensus.csv"))
tblimgt=as.matrix(read.delim("./IMGT/4_IMGT-gapped-AA-sequences.txt",check.names=F))
    idIMGTs=sub("\\_consensus.*","",tblimgt[,2])

## D/J/V primer sequences
d1="GCTTATCTGGTGGTTTCTTCCAGC"
d2="GTAGGCACCTGTGGGGAAGAAACT"
j1="GGTATAAGGTAAGACAGAGTCGTCCCT"
j2="TTAGGGAATCTCCGGGAGGGAAA"
dv34="AGAGTCGGTGGTGCAACTGAACCT"      # TRBV31
## germDJ1, from D1 primer to Trbj1-1
dj1="GCTTATCTGGTGGTTTCTTCCAGCCCTCAAGGGGTAGACCTATGGGAGGGTCCTTTTTTGTATAAAGCTGTAACATTGTGGGGACAGGGGGCCACGGTGATTCAATTCTATGGGAAGCCTTTACAAAAACCATTCTGTCTGTCCCAAGGCCCACAGTCCTCTACTCTTCTCCTGGGTCGCAAGCCTAGATGTGTTTAGATCCAGAATGCTTTCACGTCACCTCTGCAGCCTGCTAGGCCAAGATTGTGGAGAAGAATGGGGCCTTACCAGAGCAGGAGCCTCCTACACTGAATGAACATAAGAAGAATGGCCTAGTGGCCCTAGCAGCAAAGGATGCAGGGAAATCCTCTTCTGCCAACGTCCATAAAACAGACCACCATCAGTGGATAGGTGAGCCAGAGGTGAACAGTTCTACAGGTTTGTTTTGTTTTAAAAGATCTCAGCTCTTGATGAATATCATCATAGGACCAGGAAAACTGTTCTCTGACAGGCTACCTCACTTTGATGGACCTTGAGAGAGGAACAATATATATGGGGGATTCTTCACAAAAGGGATGTAAGTTCCACTGGAGGGAATCTACCATGTTTGACATTGCCACAAGTCTATCGTCAAAATGTTCTTTTGTGCAATACCTGTGATGCACACAAAGCGATTACTCCTCCTATGGTCCTCTAGATACTAGCTGGGAAGAGCCTCTCTTCACCCCTTAAGATTATTTTTCTCCTCATCCTATGGCACTGTGCAAACACAGAAGTCTTCTTTGGTAAAGGAACCAGACTCACAGTTGTAG"
## germDJ2, from D2 primer to Trbj2-1
dj2="GTAGGCACCTGTGGGGAAGAAACTTTTTTGTATCACGATGTAACATTGTGGGGACTGGGGGGGCCACAATGATTCAACTGGAAGAGGTGCTTTTACAAAAAGCTCTACCCAAAAAAACACCAACTCTCAGCCCAAAAAATGCATTTAACATGTGAGGAGGAGTCTATGTGAGTGGACTCACAAGGTCTTATAACATCTATGCATCTTCTTGCCCTAGCAAGTTTCCCACGAGTGTGAAAGGCTATCCAGGGAAGTGGGACCACCATTCAGCCCTTATCATGTTAACCATCAGCAAAGTATCCGTACATGCATGGGAAATTATGAACTTATTAAAAAGCATTAAGTAGAAAAATAAACCCTTTGTGGGAAAGAAGTTACAAGTCCTTTCTCATTGGAGGCTACCTGGACACCAACAGGAAGCCTTCAGATTTCAAATTCACTTTGTCAAACCCTGAGTGAATAGATGGATATCCGTTCCCAAGCCAAAAGTGGTATCTCCTCCTCAATTTGAGATCGGCCTCATGCAAGGTCAAGATTGCCTTACCAGTTCTGGAGGTAGATGGAGAATGTGAGTAACCCCGGGTCTGTATTGAGGAAGGTGAGGAAAGAGGAAGAATTCTTGGTAGCCCTTTTCTGCTGTGTAACTATGCTGAGCAGTTCTTCGGACCAGGGACACGACTCACCGTCCTAG"

## preparation for Code to Tag ##################################################################
ListCode=c("1011000__", "0100110__", "0011000__", "0010010__", "0000110__", "0001100__", "0001000unproductive", "0000010unproductive", "0001001unproductive", "0000011unproductive", "0001000productive", "0000010productive", "0001001productive", "0000011productive")
ListTag=c("D1_germ_J1", "D2_germ_J2", "D1_x_J1", "D1_x_J2", "D2_x_J2", "D2_x_J1", "uV_x_J1", "uV_x_J2", "uV31_x_J1", "uV31_x_J2", "pV_x_J1", "pV_x_J2", "pV31_x_J1", "pV31_x_J2")

tblcodetag <- matrix(ListTag, nrow=length(ListTag), ncol=1)
rownames(tblcodetag)=ListCode
write.csv(tblcodetag,file="table42_tblcodetag.csv")

#################################################################################################
## output
header=c("Cluster","Num_seqs","Con_len","min1","min2","min3","min4","min5",
	"germDJ1","germDJ2","D1","J1","D2","J2","DV34","Functionality","Vregion","Jregion","Code","Tag","DV","rmJn","CDR3aa","bp", "nt")
output=matrix("__",nrow=dim(ccs)[1],ncol=length(header))
colnames(output)=header

rownames(output)=ccs[,1]

## consensus
output[,1:8]=ccs[,1:8]

## IMGT information
for(j in which(idIMGTs %in% ccs[,1])){ 
    idIMGT=sub("\\_consensus.*", "", tblimgt[j,2])

# Functionality
	output[idIMGT,"Functionality"]=tblimgt[j,"V-DOMAIN Functionality"]
# TRBV
    trbv<-sub("\\*.*", "", tblimgt[j,"V-GENE and allele"])
	output[idIMGT,"Vregion"]=sub("Musmus ", "", trbv)
# TRBJ
    trbj=sub("\\*.*", "", tblimgt[j,"J-GENE and allele"])
	output[idIMGT,"Jregion"]=sub("Musmus ", "", trbj)
# CDR3 aa
    output[idIMGT,"CDR3aa"]=(paste0("...", tblimgt[j,"CDR3-IMGT"], "..."))
}

## a small function to count the number of matches
count_grep <- function(pattern, query){
	tmp=unlist(gregexpr(pattern,query,fixed=T)) 
	if(any(tmp<0)){ # fixed=T; otherwise long pattern leads to error...omg R..
		output=0
	}else{
		output=length(tmp)
	}
	return(output)
}

## a small function to get last 6 characters
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

## presence of germ/primer sequence
for(i in 1:dim(output)[1]){
	fafile=paste0("fa_consensus/",output[i,1],"_consensus.fa")
    fa=readLines(fafile)[2]
# code
    output[i,"germDJ1"]=count_grep(dj1,fa)
	output[i,"germDJ2"]=count_grep(dj2,fa)
	output[i,"D1"]=count_grep(d1,fa)
	output[i,"J1"]=count_grep(j1,fa)
	output[i,"D2"]=count_grep(d2,fa)
	output[i,"J2"]=count_grep(j2,fa)
	output[i,"DV34"]=count_grep(dv34,fa)
    code=paste0(output[i,9], output[i,10], output[i,11], output[i,12], output[i,13], output[i,14], output[i,15], output[i,16], collapse=" ")
	output[i,"Code"]=code
    
    for(j in which(code %in% ListCode)){
        output[i,"Tag"]=tblcodetag[code,1]
    }

# DV
    dv1=sub("_cluster.*", "__________", output[i,1])
    dv2=sub("DV", "__________DV", dv1)
    dv3=substr(dv2, 16, 31)
    dv4=gsub("_", "", dv3)
    dv5=as.numeric(sub("DV", "", dv4))
    output[i,"DV"]=dv5
# rmJn, sequence
    if (dv5 >= 3) {
        rmJnfile=paste0("fa_consensus_rmJn_VDJ/",output[i,1],"_consensus.rmJn.fa")
    } else {
        rmJnfile=paste0("fa_consensus_rmJn_DJ/",output[i,1],"_consensus.rmJn.fa")
    }
    rmJnName=readLines(rmJnfile)[1]
    rmJnseq=readLines(rmJnfile)[2]
    output[i,"rmJn"]=substrRight(rmJnName, 6)
    output[i,"nt"]=tolower(rmJnseq)    
    output[i,"bp"]=nchar(rmJnseq)    
}

row.names(output)=NULL
write.csv(output,file="table42_Cluster.csv")

