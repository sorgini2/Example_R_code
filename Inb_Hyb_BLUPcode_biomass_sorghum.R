#Inb_Hyb_BLUP for biomass sorghum
rm(list=ls())
setwd("/Volumes/themes/EBI/Database/")
pheno=read.csv("biomass_sorghum_phenotypes.csv",header=T,stringsAsFactors=F)
str(pheno)

#fix missing, seg values in pheno$mr
table(pheno$mr)
#pheno[which(pheno$mr==""),7]=NA
pheno[which(pheno$mr=="2,4"),7]=NA
pheno$mr=as.numeric(pheno$mr)

#change 12EF moisture values to moist new values
pheno[pheno$location=="12EF",24]=pheno[pheno$location=="12EF",25]
pheno_12EF=pheno[pheno$location=="12EF",]
pheno_13EF=pheno[pheno$location=="13EF",]


####look at spatial/nitrogen effects on 13EF checks
checks13=pheno_13EF[pheno_13EF$name=="Pacesetter",]
lm1=lm(checks13[,26]~checks13[,27])
plot(checks13[,27],checks13[,26])
abline(lm1)
summary(lm1) #r2 =0.38 (nitrogen explains 38% of check variance in yield)
normalized_yield=(checks13[,26]-mean(checks13[,26],na.rm=T))/sd(checks13[,26],na.rm=T) #substract mean, divide by stdev
colvec=ifelse(normalized_yield>0,"lightblue","orange")
cexvec=abs(normalized_yield)*8
checks13$row=as.numeric(gsub("_.+","",checks13$row)) #change row format to numeric
labelvec=ifelse(checks13[,27]==1,"G",ifelse(checks13[,27]==1.5,"GY",ifelse(checks13[,27]==2,"Y",NA)))
pdf("pdfs/Pacesetter_checks_13EF_yield_nitrogen_new.pdf",width=9,height=7)
plot(checks13$range,checks13$row,col=colvec,cex=cexvec,pch=19,cex.lab=1.2,cex.axis=1.2,ylab="Row",xlab="Range")
text(checks13$range,checks13$row,labels=labelvec,cex=2)

abline(v=5,lty=3)
abline(v=9,lty=3)
abline(v=13,lty=3)
abline(v=17,lty=3)
abline(v=21,lty=3)
abline(v=25,lty=3)
abline(v=29,lty=3)
abline(v=33,lty=3)
abline(v=37,lty=3)
abline(h=24,lty=3)
dev.off()

abline(h=16,lty=3)dev.off()

#select columns to include
phenos=apply(pheno[,c(7,11:20,24,26)],2,as.numeric)
phenos12=apply(pheno_12EF[,c(7,11:20,24,26)],2,as.numeric)
phenos13=apply(pheno_13EF[,c(7,11:20,24,26)],2,as.numeric)

GENO = as.factor(pheno$name)
BLOCK = as.factor(pheno$Block)
LOC = as.factor(pheno$location)

GENO12 = as.factor(pheno_12EF$name)
BLOCK12 = as.factor(pheno_12EF$Block)


GENO13 = as.factor(pheno_13EF$name)
BLOCK13 = as.factor(pheno_13EF$Block)


entries=names(table(GENO))
entries12=names(table(GENO12))
entries13=names(table(GENO13))

blup_matrix=matrix(NA,length(table(GENO)),ncol(phenos))
blup_matrix12=matrix(NA,length(table(GENO12)),ncol(phenos12))
blup_matrix13=matrix(NA,length(table(GENO13)),ncol(phenos13))

library(lme4)

for (i in 1:ncol(phenos)){
  temppheno=phenos[,i]
  temppheno12=phenos12[,i]
  temppheno13=phenos13[,i]
  model = lmer(temppheno ~ (1|LOC) + (1|LOC:BLOCK) + (1|GENO)) 
  model12 = lmer(temppheno12 ~ (1|BLOCK12) + (1|GENO12))
  model13 = lmer(temppheno13 ~ (1|BLOCK13) + (1|GENO13))
  blup=ranef(model)
  blup12=ranef(model12)
  blup13=ranef(model13)
  #blup_matrix[,i]=round(blup$GENO[,1],4) #doesn't work when there are NAs for some entries
  for (j in 1:length(entries)){
    if (entries[j]%in%rownames(blup$GENO)==F) blup_matrix[j,i]=NA else{
      blup_matrix[j,i]=blup$GENO[which(rownames(blup$GENO)==entries[j]),1]}
  }
  for (j in 1:length(entries12)){
    if (entries12[j]%in%rownames(blup12$GENO12)==F) blup_matrix12[j,i]=NA else{
      blup_matrix12[j,i]=blup12$GENO12[which(rownames(blup12$GENO12)==entries12[j]),1]}
  }
  for (j in 1:length(entries13)){
    if (entries13[j]%in%rownames(blup13$GENO13)==F) blup_matrix13[j,i]=NA else{
      blup_matrix13[j,i]=blup13$GENO13[which(rownames(blup13$GENO13)==entries13[j]),1]}
  }
  rawmean=tapply(temppheno,GENO, na.rm=T, mean)
  rawmean12=tapply(temppheno12,GENO12, na.rm=T, mean)
  rawmean13=tapply(temppheno13,GENO13, na.rm=T, mean)
  pdf(paste("BLUPs/12EF13EF_",colnames(phenos)[i],"_BLUPs_vs_Means.pdf",sep=""))
  plot(rawmean, blup_matrix[,i], xlab="Mean Data", ylab="BLUPs", col= "blue",main=colnames(phenos)[i])
  dev.off()
  pdf(paste("BLUPs/12EF_",colnames(phenos)[i],"_BLUPs_vs_Means.pdf",sep=""))
  plot(rawmean12, blup_matrix12[,i], xlab="Mean Data", ylab="BLUPs", col= "blue",main=colnames(phenos)[i])
  dev.off()
  pdf(paste("BLUPs/13EF_",colnames(phenos)[i],"_BLUPs_vs_Means.pdf",sep=""))
  plot(rawmean13, blup_matrix13[,i], xlab="Mean Data", ylab="BLUPs", col= "blue",main=colnames(phenos)[i])
  dev.off()
  }

blup_matrix=data.frame(names(table(GENO)),blup_matrix)
colnames(blup_matrix)=c("Taxa",colnames(phenos))
write.table(blup_matrix,file="BLUPs/All_Inbred_BLUPs_new.csv",quote=F,sep=",",row.names=F)

blup_matrix12=data.frame(names(table(GENO12)),blup_matrix12)
colnames(blup_matrix12)=c("Taxa",colnames(phenos12))
write.table(blup_matrix12,file="BLUPs/12EF_Inbred_BLUPs_new.csv",quote=F,sep=",",row.names=F)

blup_matrix13=data.frame(names(table(GENO13)),blup_matrix13)
colnames(blup_matrix13)=c("Taxa",colnames(phenos13))
write.table(blup_matrix13,file="BLUPs/13EF_Inbred_BLUPs_new.csv",quote=F,sep=",",row.names=F)


