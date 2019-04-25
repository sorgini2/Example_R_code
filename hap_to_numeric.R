setwd("C:/Users/sorgini2/Desktop/hap_to_num")
directory="C:/Users/sorgini2/Desktop/hap_to_num/"
file_list1=list.files(directory)
i=1
data=read.table(paste(directory,file_list1[i],sep=""),header=T,sep="\t",comment.char="",stringsAsFactors=F) # lines containing "#" ignored unless you turn off comment.char

name="ALLZEA_chr1_2_"
het_thresh=0.2
maf_thresh=0

alldata=matrix(NA,0,ncol(data)-9)
for (i in 1:length(file_list1)){
  data=read.table(paste(directory,file_list1[i],sep=""),header=T,sep="\t",comment.char="",stringsAsFactors=F) 
  #if (nrow(data)==1) next
  indels=grepl("-",data[,2])
  data=data[indels==F&nchar(data[,2])>1,]
  major=substr(data[,2],1,1)
  minor=substr(data[,2],3,3)
  if("REFERENCE" %in% colnames(data)==T) {ref=data[,which(grepl("REFERENCE",colnames(data))==T)]; alt=ifelse(ref==major,minor,major)} #ref=0,alt=1
  if("REFERENCE" %in% colnames(data)==F) {ref=major; alt=minor} #major=0,minor=1
  snps=as.matrix(data[,12:ncol(data)])
  for (j in 1:nrow(snps)){
    snps[j,]=gsub(ref[j],0,snps[j,])  
    snps[j,]=gsub(alt[j],2,snps[j,])
    snps[j,][snps[j,]%in%c("R","Y","S","W","K","M")]=1
    snps[j,][snps[j,]=="N"]=NA
    snps[j,][snps[j,]=="X"]=NA
  }
  snps=cbind(data[,3:4],snps)
  snps=as.matrix(snps)
  alldata=rbind(alldata,snps)
  print(i)
}

alldata=apply(alldata,2,as.numeric)
snp_locs=alldata[,1:2]
alldata=alldata[,3:ncol(alldata)]
taxa=sapply(strsplit(colnames(alldata),"\\."),"[[",1) # useful piece of code for extracting taxa names
snps012=alldata; rm(alldata)
snp_info=snp_locs; rm(snp_locs)
colnames(snps012)=taxa

#apply filters
nas=apply(snps012,1,function(x){sum(is.na(x))})
snps012=snps012[nas==0,]; snp_info=snp_info[nas==0,]

if (het_thresh>0){
  hets=apply(snps012,1,function(x){sum(x==1,na.rm=T)/ncol(snps012)})
  hist(hets)
  snps012=snps012[hets< het_thresh,]; snp_info=snp_info[hets< het_thresh,]
}
if (maf_thresh>0){
  maf=apply(snps012,1,function(x){mean(x,na.rm=T)/2})
  snps012=snps012[maf > maf_thresh & maf < (1-maf_thresh),]; snp_info=snp_info[maf > maf_thresh & maf < (1-maf_thresh),]
}
setwd(directory)
save(snps012,file=paste(name,"_snps012.RData",sep=""))
save(snp_info,file=paste(name,"_info.RData",sep=""))
save(taxa,file=paste(name,"_taxa.RData",sep=""))

