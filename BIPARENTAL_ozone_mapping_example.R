## RESPONSE BIPARENTAL QTL MAPPING
## STEPWISE REGRESSION
## 
##NILs WITH B73 AND MO17 RUN SEPARATE 
#test models and extract geno_effects
data$ring=as.factor(data$ring)
ele=data[data$treatment=="elevated",]
amb=data[data$treatment=="ambient",]
lm1=lm(ele$value~ele$ring+ele$ringset%in%ele$ring+ele$genotype)
lm2=lm(ele$value~ele$ring+ele$ringset+ele$genotype)
AIC(lm1)#1231.037
AIC(lm2)#1272.1
lm3=lm(amb$value~amb$genotype)
geno_effects_amb=lm3$coefficients[grepl("genotype",names(lm3$coefficients))]
geno_effects_ele=lm1$coefficients[grepl("genotype",names(lm1$coefficients))]
taxa=gsub("ele.genotype","",names(geno_effects_ele))

out=data.frame(taxa,geno_effects_amb,geno_effects_ele)
plot(out[,2],out[,3])
#out[out[,1]=="Mo17",1]="mo17"

mo_out=out[grepl("^m",out[,1])==T|grepl("^M",out[,1])==T,]
response=round(mo_out[,3]-mo_out[,2],2)
mo_out=data.frame(mo_out,response)
b73_out=out[grepl("^b",out[,1])==T|grepl("^B",out[,1])==T,]
response2=round(b73_out[,3]-b73_out[,2],2)
b73_out=data.frame(b73_out,response2)

#import genotypes/map
geno=read.csv("b73_mo17_nils.csv")

##Scan for QTL#1

geno2=geno[,colnames(geno)%in%mo_out[,1]]
geno3=geno[,colnames(geno)%in%b73_out[,1]]
Mo17=rep(2,nrow(geno2))
B73=rep(0,nrow(geno3))
geno2=cbind(geno2,Mo17)
geno3=cbind(geno3,B73)
colnames(geno2)==mo_out[,1]
colnames(geno3)==b73_out[,1]
geno2=as.matrix(geno2)
geno3=as.matrix(geno3)
geno2=apply(geno2,2,as.numeric)
geno3=apply(geno3,2,as.numeric)
pvals=c()
pvals2=c()
for (i in 1:nrow(geno2)){
  pvals[i]= anova(lm(mo_out$response~geno2[i,]))[1,5]
  pvals2[i]= anova(lm(b73_out$response~geno3[i,]))[1,5]
  print(i)
}

colvec=ifelse(geno[,1]%in%c(1,3,5,7,9),"orange","blue")
#mo17 Results
#pdf("mo17_nils_results.pdf")
plot(-log10(pvals),col=colvec,pch=19, main="Maize [O3] 2016 Leaf Damage", font=2)
#plot(-log10(pvals)[1500:1900],col=colvec,pch=19)
plot(-log10(pvals),col=colvec,pch=19, main="",xlab="", font=2, xaxt="n", ylab=expression(bold("-log10(pvals)")))

abline(h=6.4,lwd=3, col="black", lty=1)
#plot(-log10(pvals2),col=colvec,pch=19)
#dev.off()
geno[which.min(pvals),] #chr2 160.938 Mb
geno[1782,]
anova(lm(mo_out$response~geno2[1782,]))

-log10(2.467e-08) # peak pvalue

pvals_step2=c()
for (i in 1:nrow(geno2)){
  pvals_step2[i]= anova(lm(mo_out$response~geno2[1782,]+geno2[i,]))[2,5]
  print(i)
}

geno[which.min(pvals_step2),] #4985   6 160297015 
min(pvals_step2, na.rm=TRUE)
#[1] 0.004546643
geno[4985,]
plot(-log10(pvals_step2),col=colvec,pch=19)#not over significance threshold

##step3
pvals_step3=c()
for (i in 1:nrow(geno2)){
  pvals_step3[i]= anova(lm(mo_out$response~geno2[1782,]+geno2[4985,]+geno2[i,]))[3,5]
  print(i)
}
geno[which.min(pvals_step3),] #  2 198871298.1 161.18
min(pvals_step3, na.rm=TRUE)
#[1] 0.002816488
geno[1896,]
plot(-log10(pvals_step3),col=colvec,pch=19)#not over significance threshold

### WITH RECURRENT PARENT AS COVARIATE

ele=data[data$treatment=="elevated",]
amb=data[data$treatment=="ambient",]
#lm1=lm(ele$value~ele$ring+ele$ringset%in%ele$ring+ele$genotype)
lm1=lm(ele$value~ ele$year + ele$ringset%in%ele$year + ele$entryset%in%ele$year + ele$genotype)
#lm2=lm(ele$value~ele$ring+ele$ringset+ele$genotype)
lm2=lm(ele$value~ ele$year + ele$ringset + ele$genotype)
AIC(lm1)#2255.975
AIC(lm2)#2252.125
lm3=lm(amb$value~amb$genotype)
geno_effects_amb=lm3$coefficients[grepl("genotype",names(lm3$coefficients))]
geno_effects_ele=lm2$coefficients[grepl("genotype",names(lm2$coefficients))]
taxa=gsub("ele.genotype","",names(geno_effects_ele))

out=data.frame(taxa,geno_effects_amb,geno_effects_ele)
plot(out[,2],out[,3])
#out[out[,1]=="Mo17",1]="mo17"

response=round(out[,3]-out[,2],2)
out=data.frame(out,response)
write.csv(out, file = "out.csv")
###ADD RP to file; and add extra column place holder
out=read.csv("test.csv")
out$RP=as.factor(out$RP)
str(out$RP)
out$taxa=as.character(out$taxa)
str(out$taxa)
out=data.frame(out)
#import genotypes/map
geno=read.csv("b73_mo17_nils.csv")
#View(geno)
##Scan for QTL#1 RESPONSE
#geno2=geno[,colnames(geno)%in%out[,2]]
geno2=geno[,colnames(geno)%in%out[,3]]
#geno3=geno[,colnames(geno)%in%out[,2]]
Mo17=rep(2,nrow(geno2))
B73=rep(0,nrow(geno2))
geno2=cbind(geno2,Mo17)
geno2=cbind(geno2,B73)
#geno3=cbind(geno3,B73)

colnames(geno2)==out[,3]
#colnames(geno3)==out[,2]
geno2=as.matrix(geno2)
#geno3=as.matrix(geno3)
geno2=apply(geno2,2,as.numeric)#50 warnings produced?
#geno3=apply(geno3,2,as.numeric)
pvals=c()
#pvals2=c()
for (i in 1:nrow(geno2)){
  pvals[i]= anova(lm(out$response~out$RP+geno2[i,]))[2,5]#column 2 becasue 1 is covariate
  #pvals2[i]= anova(lm(out$response~out$RP+geno3[i,]))[1,5]
  print(i)
}
colvec=ifelse(geno[,1]%in%c(1,3,5,7,9),"orange","blue")
plot(-log10(pvals),col=colvec,pch=19, main="", font=2)
min(pvals,na.rm=TRUE)
#[1] 1.00333e-13
geno[which.min(pvals),] #1782   2 160938561.9
geno[1782,]
anova(lm(out$response~out$RP+geno2[1782,]))
####
pvals_step2=c()
for (i in 1:nrow(geno2)){
  pvals_step2[i]= anova(lm(out$response~out$RP+geno2[1782,]+geno2[i,]))[3,5]
  print(i)
}
geno[which.min(pvals_step2),] #5519   7 150513706 137.58 
min(pvals_step2, na.rm=TRUE)
#[1] 0.0007323286
geno[5519,]
plot(-log10(pvals_step2),col=colvec,pch=19)#not over significance threshold

##step3
pvals_step3=c()
for (i in 1:nrow(geno2)){
  pvals_step3[i]= anova(lm(out$response~out$RP+geno2[1782,]+geno2[5519,]+geno2[i,]))[4,5]
  print(i)
}
geno[which.min(pvals_step3),] #876   1 218987992.5 249.48 
min(pvals_step3, na.rm=TRUE)
#  0.001477056
geno[876,]
plot(-log10(pvals_step3),col=colvec,pch=19)#not over significance threshold

###BY ENVIRONMENT (ele)
pvals=c()
for (i in 1:nrow(geno2)){
  pvals[i]= anova(lm(out$geno_effects_ele~out$RP+geno2[i,]))[2,5]
  print(i)
}

colvec=ifelse(geno[,1]%in%c(1,3,5,7,9),"orange","blue")
#Mo17 results ELEVATED
plot(-log10(pvals),col=colvec,pch=19, main="")
geno[which.min(pvals),] ##chr2 160.938 Mb
geno[1782,]
anova(lm(out$geno_effects_ele~out$RP+geno2[1782,]))
boxplot(out$geno_effects_ele~out$RP+geno2[1782,],ylab="",xlab="# of Alleles", main="mo_out$geno_effects_ele~geno2[step1]", cex.lab=1,cex.axis=1, las=1)

##step2 BY ENV ele
pvals_step2=c()
for (i in 1:nrow(geno2)){
  pvals_step2[i]= anova(lm(out$geno_effects_ele~out$RP+geno2[1782,]+geno2[i,]))[3,5]
  print(i)
}

geno[which.min(pvals_step2),] #5519   7 150513706 137.58
min(pvals_step2, na.rm=TRUE)
# 0.001540661
geno[5519,]
plot(-log10(pvals_step2),col=colvec,pch=19)#not over significance threshold of 5.3
##step 3
pvals_step3=c()
for (i in 1:nrow(geno2)){
  pvals_step3[i]= anova(lm(out$geno_effects_ele~out$RP+geno2[1782,]+geno2[5519,]+geno2[i,]))[4,5]
  print(i)
}
geno[which.min(pvals_step3),] #876   1 218987992.5 249.48  
min(pvals_step3, na.rm=TRUE)
# 0.001705912
geno[876,]
plot(-log10(pvals_step3),col=colvec,pch=19)#not over


