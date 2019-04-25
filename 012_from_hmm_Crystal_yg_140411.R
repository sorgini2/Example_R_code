#hmm parent format to 012
#Maps genomic region responsible for YG3 mutant phenotype
  setwd("~/illum/maize/hmm2")
  directory="~/illum/maize/hmm2/"
    file_list1=list.files(directory)
    i=1
    data=read.table(paste(directory,file_list1[i],sep=""),header=T,sep="\t",comment.char="",stringsAsFactors=F) # lines containing "#" ignored unless you turn off comment.char
    alldata=matrix(NA,0,ncol(data)-9)
    for (i in 1:length(file_list1)){
      data=read.table(paste(directory,file_list1[i],sep=""),header=T,sep="\t",comment.char="",stringsAsFactors=F)  
      snps=as.matrix(data[,12:ncol(data)])
      snps[snps=="A"]=0
      snps[snps=="M"]=1
      snps[snps=="C"]=2
      snps[snps=="N"]=NA
      snps=apply(snps,2,as.numeric)
      snps=cbind(data[,3:4],snps)
      snps=apply(snps,2,as.numeric)
      alldata=rbind(alldata,snps)
      print(i)
    }
  alldata=apply(alldata,2,as.numeric)
  alldata=alldata[order(alldata[,1],alldata[,2]),]
  snp_locs=alldata[,1:2]
  alldata=alldata[,3:ncol(alldata)]
  taxa=sapply(strsplit(colnames(alldata),"\\."),"[[",1)
  

snps012=alldata; rm(alldata); rm(snps); rm(data)
colnames(snps012)=taxa

pheno=ifelse(grepl("Gree",taxa),0,ifelse(grepl("B73",taxa),0,1))
pheno[177]=0 #REFERENCE GENOME
table(snps012[,1])

pvals=c()
for (i in 1:nrow(snps012)){
  pvals[i]=anova(lm(pheno~snps012[i,]))[1,5]
  print(i)
}

  input=snp_locs
  chroms=names(table(input[,1]))
  chrom_lengths=c()
  for (i in 1:length(chroms)){
    chrom_lengths[i]=max(input[input[,1]==chroms[i],2])
  }
  
  chrom_cumul=c()
  for (i in 1:length(chrom_lengths)){
    chrom_cumul[i]=sum(chrom_lengths[1:i])
  }
  chrom_cumul=c(0,chrom_cumul)
  
  cumul=c()
  for (i in 1:nrow(input)){
    cumul[i]=chrom_cumul[input[i,1]]+input[i,2]
  }
  
  
colvec=ifelse(snp_locs[,1]%in%c(1,3,5,7,9),"black","grey")
pdf("yg_whole_genome.pdf",width=10,height=6)
plot(cumul,-log10(pvals),xaxt="n",xlab="",col=colvec)
dev.off()
pdf("yg_closeup.pdf",width=6,height=4)
plot(snp_locs[2420:2460,2]/1000000,-log10(pvals[2420:2460]),type="l",xlab="Mb on chrom 5",ylab="-log10(pvals")
dev.off()
snp_locs[which.min(pvals),] #2437

out1=cbind(snp_locs[2420:2460,],pvals[2420:2460])
out2=snps012[2420:2460,]
out3=data.frame(out1,out2)
setwd("/Volumes/themes/EBI/pjb_lab/Database/Crystal")
write.csv(out3,file="yg_interval_c5_140429.csv",quote=F,row.names=F)
