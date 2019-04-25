#calculate basic stats (heterozygosity, missing)
#for duplicated samples, calculate concordance
#for low-concordance samples, use SC-PRE kinship to find bad samples
#combine duplicated samples
rm(list=ls(all=TRUE))
load("/Volumes/themes/EBI/pjb_lab/Database/R/illum/SORGSNPS/snps012.RData")

genos=snps012[,7:ncol(snps012)]
snpstats=snps012[,1:6]
rm(snps012)
#split colnames into sample, plate, and well names
id=sapply(strsplit(colnames(genos),"\\."),"[[",1)
lane=paste(sapply(strsplit(colnames(genos),"\\."),"[[",2),sapply(strsplit(colnames(genos),"\\."),"[[",3),sep="_")
lib=sapply(strsplit(colnames(genos),"\\."),"[[",4)
well=sapply(strsplit(colnames(genos),"\\."),"[[",5)

stats=cbind(id,lane,lib,well)
genos=genos[,grepl("SC",stats[,1])]
stats=stats[grepl("SC",stats[,1]),]


missing_by_sample=apply(genos,2,function(x){sum(is.na(x))})
prop_hets_by_sample=c()
for (i in 1:ncol(genos)){
  temp=genos[,i][!is.na(genos[,i])]
  prop_hets_by_sample[i]=sum(temp==1)/sum(table(temp))
}
#avg pairwise concordance?
id=stats[,1]
dups=names(table(id)[table(id)>1])
number_dups=table(id)[table(id)>1]
name_labels=paste(dups,"(",number_dups,")",sep="")

nonconc_rate=c()
libs_involved=c()
wells_involved=c()
for (i in 1:length(dups)){
  temp=genos[,id==dups[i]]
  libs_involved[i]=paste(lib[id==dups[i]],collapse=",")
  wells_involved[i]=paste(well[id==dups[i]],collapse=",")
  nonconc_mat=matrix(NA,ncol(temp),ncol(temp))
  for (j in 1:(ncol(temp)-1)){
    for (k in (j+1):ncol(temp)){
      temp_vec1=temp[,j][!is.na(temp[,j])&!is.na(temp[,k])]
      temp_vec2=temp[,k][!is.na(temp[,j])&!is.na(temp[,k])]
      nonconc_mat[j,k]=sum(ifelse(temp_vec1!=temp_vec2,1,0))/length(temp_vec1)
    }
  }
  nonconc_rate[i]=mean(nonconc_mat[!is.na(nonconc_mat)])
  }

pdf("/Volumes/themes/EBI/pjb_lab/Database/R/illum/SORGSNPS/sc_lines_nonconc_rate.pdf")
temp1=cbind(name_labels,libs_involved,wells_involved,nonconc_rate)
temp1=temp1[order(temp1[,4]),]
boxplot((missing_by_sample/nrow(genos))*100~lib,names=1:9,xlab="library",ylab="% missing data by sample",main="Missing data across libraries")
boxplot(prop_hets_by_sample~lib,names=1:9,xlab="library",ylab="% heterozygous by sample",main="Heterozygosity across libraries")
par(las=2)
barplot(as.numeric(temp1[,4]),names=temp1[,1],cex.names=0.45,ylab="Non-concordance rate",xlab="Line (wells)",main="Mean pairwise non-concordance rate of inbreds \n genotyped in multiple lanes/wells")
dev.off()

#decide whether to combine samples based on their concordance rate
good_match=dups[nonconc_rate<0.07]
bad_match=dups[nonconc_rate>0.07]
unchanged=genos[,id%in%dups==F] #unduplicated samples stay the same
colnames(unchanged)=id[id%in%dups==F]


#merge duplicated samples with low nonconc rate
merged=matrix(NA,nrow(genos),length(good_match))
colnames(merged)=good_match
for (a in 1:length(good_match)){
  temp=genos[,id==good_match[a]]
  test=apply(temp,1,function(x){round(mean(x[!is.na(x)]),0)})
  merged[,a]=as.numeric(gsub(NaN,NA,test))
}

sc_genos=cbind(unchanged,merged)
sc_genos=sc_genos[,order(colnames(sc_genos))]

#filter snps by %missing and maf
missing_by_snp=apply(sc_genos,1,function(x){sum(is.na(x))})
percent_missing=missing_by_snp/ncol(sc_genos)
maf_temp=apply(sc_genos,1,function(x){sum(x[!is.na(x)])})
maf=maf_temp/((ncol(sc_genos)*2)-missing_by_snp)

table(percent_missing<0.8&maf>0.1) #9632 SNPs

sc_genos=sc_genos[percent_missing<0.8&maf>0.1,]
sc_snpstats=snpstats[percent_missing<0.8&maf>0.1,]
sc_genos=as.matrix(sc_genos)
sc_genos[sc_genos==1]=0.5
sc_genos[sc_genos==2]=1
save(sc_genos,file="/Volumes/themes/EBI/pjb_lab/Database/R/illum/SORGSNPS/sc_genos.RData")
save(sc_snpstats,file="/Volumes/themes/EBI/pjb_lab/Database/R/illum/SORGSNPS/sc_snpstats.RData")

##change snps (not-imputed) into beagle format
load("illum/SORGSNPS/sc_genos.RData")
load("illum/SORGSNPS/sc_snpstats.RData")



col1=rep("M",nrow(sc_genos_new))
snp_ids=(paste(sc_snpstats[,1],sc_snpstats[,2],sep="_"))

sc_genos_beagle=matrix(NA,nrow(sc_genos),ncol(sc_genos)*2)
colnames(sc_genos_beagle)=rep(colnames(sc_genos),2)
for (i in 1:ncol(sc_genos)){
  sc_genos_beagle[,i]=sc_genos[,i]
  sc_genos_beagle[,i][sc_genos_beagle[,i]==0]="A"
  sc_genos_beagle[,i][sc_genos_beagle[,i]==0.5]="A"
  sc_genos_beagle[,i][sc_genos_beagle[,i]==1]="G"
  sc_genos_beagle[,i+ncol(sc_genos)]=sc_genos[,i]
  sc_genos_beagle[,i+ncol(sc_genos)][sc_genos_beagle[,i+ncol(sc_genos)]==0]="A"
  sc_genos_beagle[,i+ncol(sc_genos)][sc_genos_beagle[,i+ncol(sc_genos)]==0.5]="G"
  sc_genos_beagle[,i+ncol(sc_genos)][sc_genos_beagle[,i+ncol(sc_genos)]==1]="G"
  print(i)
  }
#test=cbind(sc_genos[,i],sc_genos_beagle[,i], sc_genos_beagle[,i+ncol(sc_genos)])
sc_genos_beagle=sc_genos_beagle[,order(colnames(sc_genos_beagle))]

sc_genos_beagle[is.na(sc_genos_beagle)]="?"
sc_genos_beagle=cbind(col1,snp_ids,sc_genos_beagle)
colnames(sc_genos_beagle)[1]="I"
write.table(sc_genos_beagle,file="/Volumes/themes/EBI/pjb_lab/Database/R/illum/SORGSNPS/sc_genos_beagle.txt",quote=F,row.names=F,sep="\t")

#source("/Volumes/themes/EBI/pjb_lab/Database/R/functions/emma.R")
#sc_k_matrix=emma.kinship(sc_genos,"additive","all")
#save(sc_k_matrix,file="/Volumes/themes/EBI/pjb_lab/Database/R/illum/SORGSNPS/k_matrix.RData")




results=emma.ML.LRT(all_phenos,goodsnps,k_matrix)
