
#need a reference set of markers with cM positions additively across all chrs
genotypes=read.csv(file="book2.csv")
#genotypes <- genotypes[-1,]
genotypes_positioned <- genotypes
#genotypes with positions combined so all positions are on one scale
for (i in 2:10){
  this_chr <- genotypes[genotypes[,2]==i,]
  earlier_chr <- genotypes_positioned[genotypes_positioned[,2]<i,]
  length <- earlier_chr[nrow(earlier_chr),3]
  this_chr[,3] <- this_chr[,3] + length 
  genotypes_positioned[genotypes_positioned[,2]==i,3] <- this_chr[,3]
}
#find the last xaxis data point positon
genotypes_positioned[1478,]
#marker chrom   cm
#1478  m1478    10 1434
#add a dummy marker position to extend the axis for changing the axis formating to dotted
C = c("m0000", 10, 1600)
test=rbind(genotypes_positioned,C)
genotypes_positioned=test

#write.csv(genotypes_positioned, file = "DataGenotypes_combinedPositions.csv")
par(mar=c(0, 0, 0, 1), mfrow=c(12,1),
    oma = c(4, 4, 0.2, 0.2))
#plot(test[1:1600,3],test[1:1600,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)

#BORON (B11)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("B11")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
abline(h=10, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(955,5,col = "black",pch="-", cex=4.5)#cml103=black #jitter position from 940 to 955
points(1174,5,col = "orange",pch="-", cex=4.5)#cml333=orange
points(1113,4,col = "grey65",pch="-", cex=4.5)#Nc358 =grey65
points(1114,8,col = "grey65",pch="-", cex=4.5)#Nc358 =grey65 #jitter senconf step QTL to 6

#Magnesium (Mg26)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("Mg26")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(1368,3,col = "black",pch="-", cex=4.5)#cml103
points(60,4,col = "black",pch="-", cex=4.5)#cml103
points(1340,5.5,col = "blue",pch="-", cex=4.5)#Tx303 #jitter position from 1326 to 1340
points(80,7,col = "blue",pch="-", cex=4.5)#Tx303 #jitter position
points(1242,5,col = "grey65",pch="-", cex=4.5)#Nc358
points(1340,8,col = "grey65",pch="-", cex=4.5)#Nc358

#Phosphorus (P31)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("P31")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(536,5,col = "black",pch="-", cex=4.5)#CMl103
points(1050,5,col = "orange",pch="-", cex=4.5)#CMl333 #jitter position from 1057 to 1050

#Sulfur (S34)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("S34")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(323,5,col = "orange",pch="-", cex=4.5)#cml333
points(873,5,col = "orange",pch="-", cex=4.5)#cml333

#Potassium (K39)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("K39")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(15,4,col = "blue",pch="-", cex=4.5)#Tx303 jitter position from 2 to 15
points(23,7,col = "grey65",pch="-", cex=4.5)#Nc358 jitter from 12 to 23

#Iron (Fe54)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("Fe54")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(852,5,col = "black",pch="-", cex=4.5)#cml103
points(15,5,col = "black",pch="-", cex=4.5)#cml103 jitter from 6 to 1
points(313,5,col = "black",pch="-", cex=4.5)#cml103

#Manganese (Mn55)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("Mn55")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(1369,5,col = "black",pch="-", cex=4.5)#cml103
points(156,5,col = "black",pch="-", cex=4.5)#cml103
points(713,5,col = "blue",pch="-", cex=4.5)#Tx303
points(1287,5,col = "grey65",pch="-", cex=4.5)#Nc358

#Cobalt (Co59)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("Co59")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(980,5,col = "blue",pch="-", cex=4.5) #Tx303
points(396,5,col = "grey65",pch="-", cex=4.5) #Nc358

#Copper (Cu63)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("Cu63")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(394,5,col = "grey65",pch="-", cex=4.5) #Nc358

#Rubidium (Rb85)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("Rb85")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(477,5,col = "grey65",pch="-", cex=4.5) #Nc358

#Strontium (Sr88)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("Sr88")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(263,5,col = "black",pch="-", cex=4.5) #cml103
points(972,5,col = "orange",pch="-", cex=4.5) #cml333

#Molybdenum (Mo98)
plot(genotypes_positioned[,3],genotypes_positioned[,2],type="n", xaxt="n",yaxt="n", xlab="",ylab="", bty='n',las=2)
axis(1, col="grey75",lty=5,labels=FALSE,lwd.tick=0)
axis(2, at=1:10, labels=c("","","","","",expression(bold("Mo98")),"","","",""), las=2, tick=FALSE)
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=0, lwd=2, col="black",lty=1)#chr1_start
abline(v=202, lwd=2, col="grey75",lty=3)#chr2_start
abline(v=357, lwd=2, col="grey75",lty=3)#chr3_start
abline(v=517, lwd=2, col="grey75",lty=3)#chr4_start
abline(v=666, lwd=2, col="grey75",lty=3)#chr5_start
abline(v=821, lwd=2, col="grey75",lty=3)#chr6_start
abline(v=931, lwd=2, col="grey75",lty=3)#chr7_start
abline(v=1069, lwd=2, col="grey75",lty=3)#chr8_start
abline(v=1197, lwd=2, col="grey75",lty=3)#chr9_start
abline(v=1322,lwd=2, col="grey75",lty=3)#chr10_start
abline(v=1434, lwd=2, col="black",lty=1)#chr10_end
abline(h=.75, lwd=2, col="black",lty=1)#chr10_end
#PLOT SIGNIFICANT QTL
points(1234,5,col = "black",pch="-", cex=4.5) #cml103

#Cadmium (Cd111)
#plot(genotypes_positioned[1:1425,3],genotypes_positioned[1:1425,2], type="n", xaxt="n",yaxt="n", xlab="",ylab="", las=2)
#axis(2, at=1:10, labels=c("","","","","",expression(bold("Cd111")),"","","",""), las=2, tick=FALSE)
#abline(v=0, lwd=2)#chr1_start
#abline(v=202, lwd=2)#chr2_start
#abline(v=357, lwd=2)#chr3_start
#abline(v=517, lwd=2)#chr4_start
#abline(v=666, lwd=2)#chr5_start
#abline(v=821, lwd=2)#chr6_start
#abline(v=931, lwd=2)#chr7_start
#abline(v=1069, lwd=2)#chr8_start
#abline(v=1197, lwd=2)#chr9_start
#abline(v=1322,lwd=2)#chr10_start
#abline(v=1434)#chr10_end
#PLOT SIGNIFICANT QTL


#LABEL CHROMSOMES
axis(1, at=102.5, labels=c(expression(bold("1"))), las=1, tick=FALSE)#chromosome1
axis(1, at=286.5, labels=c(expression(bold("2"))), las=1, tick=FALSE)#chromosome2
axis(1, at=434.5, labels=c(expression(bold("3"))), las=1, tick=FALSE)#chromosome3
axis(1, at=590.5, labels=c(expression(bold("4"))), las=1, tick=FALSE)#chromosome4
axis(1, at=745, labels=c(expression(bold("5"))), las=1, tick=FALSE)#chromosome5
axis(1, at=880, labels=c(expression(bold("6"))), las=1, tick=FALSE)#chromosome6
axis(1, at=1000, labels=c(expression(bold("7"))), las=1, tick=FALSE)#chromosome7
axis(1, at=1140, labels=c(expression(bold("8"))), las=1, tick=FALSE)#chromosome8
axis(1, at=1270, labels=c(expression(bold("9"))), las=1, tick=FALSE)#chromosome9
axis(1, at=1385, labels=c(expression(bold("10"))), las=1, tick=FALSE)#chromosome10

mtext(expression(bold("Chromosome")), side = 1, outer = TRUE, line = 2.5, cex=.7)

