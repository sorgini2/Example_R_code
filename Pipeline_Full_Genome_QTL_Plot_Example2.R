setwd("~/Dropbox/IBM_GxE_Directory")
library(plyr)
library(qtl)

#Set locations
locations <- c("NY05", "NY06", "NC06", "FL06", "MO06", "SA10", "IN09", "IN10", "NY12", "FL05")

#Set type (this is bc I had element, PC traits etc., can also just edit file name below)
type <- "El"

#Load in stepwise qtl output rdata files and crosses
stepwise_list_total <- list()
cross_total <- list()
for (z in 1:length(locations)){
  
  #change files below to your file paths
  load(file = paste("Data/RQTL_Output/", type, "/", locations[z], "_", type, "_Stepwise_List.Rdata", sep=""))
  stepwise_list_large <- stepwise_list
  #load in cross
  load(file = paste("Data/Crosses/", type, "/", type, "_RQTL_cross_", locations[z], ".RData", sep=""))
  cross_large <- cross
  cross_large <- sim.geno(cross_large)
  
  stepwise_list_total[[z]] <- stepwise_list_large
  cross_total[[z]] <- cross_large
}

#need a reference set of markers with cM positions additively across all chrs
#load in your genotypes
genotypes <- read.csv(file = "Data/Raw_Data/IBM_MAP_4000markers_forRQTL.csv", sep = ",", header = FALSE, stringsAsFactors = FALSE)
genotypes <- genotypes[-1,]
genotypes_positioned <- genotypes
#genotypes with positions combined so all positions are on one scale
for (i in 2:10){
  this_chr <- genotypes[genotypes[,2]==i,]
  earlier_chr <- genotypes_positioned[genotypes_positioned[,2]<i,]
  length <- earlier_chr[nrow(earlier_chr),3]
  this_chr[,3] <- this_chr[,3] + length 
  genotypes_positioned[genotypes_positioned[,2]==i,3] <- this_chr[,3]
}

#write.csv(genotypes_positioned, file = "Data/Genotypes_combinedPositions.csv")

load("Data/color_vector.Rdata")
#reorganize color vector or create your own - one color for each environment
color_vector <- color_vector[c(6,7,5,2,9,10,3,4,8,1)]

#for each trait go through and make data matrix
#contains trait, QTL location, and color column designating environment
#combine all the data matrices from all traits at the end
dat_full <- matrix(ncol=3)
colnames(dat_full) <- c("position", "trait_num", "env_color")
for (j in 1:nphe(cross)){
  position_env_matrix <- matrix(ncol=2)
  chr_vector <- c()
  pos_vector <- c()
  #for each location
  for (i in 1:length(locations)){
    if (length(stepwise_list_total[[i]][[j]]) != 0){
      length_before <- length(chr_vector)
      chr_vector <- c(chr_vector, stepwise_list_total[[i]][[j]]$chr)
      pos_vector <- c(pos_vector, stepwise_list_total[[i]][[j]]$pos)
      places <- c((length_before+1):(length_before + length(stepwise_list_total[[i]][[j]]$chr)))
      places <- unique(places)
      loc_matrix <- cbind(places, rep(i,length(places)))
      position_env_matrix <- rbind(position_env_matrix, loc_matrix)
    }
  }
  if (nrow(position_env_matrix)!=1){
    if (nrow(position_env_matrix)==2){
      position_env_matrix <- t(as.matrix(position_env_matrix[-1,]))
    }
    if (nrow(position_env_matrix)!=1){
      position_env_matrix <- as.matrix(position_env_matrix[-1,])
    }
    chr_vector
    pos_vector
    
    colors <- color_vector[position_env_matrix[,2]]
    qtl <- makeqtl(cross_large, chr=chr_vector, pos=pos_vector)
    
    plot(qtl, col=colors, horizontal=T, justdots=T)
    
    pos_vector_new <- pos_vector
    #find corresponding chromosome
    for (i in 1:length(pos_vector)){
      if (chr_vector[i] != 1){
        chr_vector[i]
        earlier_chr <- genotypes_positioned[as.numeric(genotypes_positioned[,2])<as.numeric(chr_vector[i]),]
        length <- earlier_chr[nrow(earlier_chr),3]
        pos_vector_new[i] <- pos_vector_new[i] + length + 250*(as.numeric(as.character(chr_vector[i]))-1)
      }
    }
    
    pos_vector
    pos_vector_new
    
    dat <- data.frame(pos_vector_new, rep(j, length(pos_vector)), colors)
    colnames(dat) <- c("position", "trait_num", "env_color")
    dat_full <- rbind(dat_full, dat)
  }
}

#the next steps jitter positions of QTL that are on top of each other
#this helps for visibility with the plot to some degree, some lines will still be kind of on top of each other
mround <- function(x,base){ 
  base*round(x/base) 
} 

dat_full_round <- dat_full
dat_full_round[,1] <- mround(dat_full[,1],5)
jittered <- jitter(dat_full_round[,1], factor=1)
jitter(dat_full_round[1:2,1])      

dat_full_round <- arrange(dat_full_round, position)
dat_full_round <- na.omit(dat_full_round)

dat_full_round_jittered <- dat_full_round
jitter_list <- list()
for (i in 1:nrow(dat_full_round)){
  if (nrow(dat_full_round[dat_full_round[,1] == dat_full_round[i,1],]) != 1 & length(unique(dat_full_round[dat_full_round[,1] == dat_full_round[i,1],][,2])) == 1){
    jitter_list[[i]] <- rownames(dat_full_round[dat_full_round[,1] == dat_full_round[i,1],])
  }
}

unique(jitter_list)

for (i in 2:length(unique(jitter_list))){ 
  x <- length(dat_full_round_jittered[unique(jitter_list)[[i]], 2])
  if (x %% 2 != 0){
    jitter_amt <- c((-(x-1)/2):((x-1)/2))
  }
  if (x %% 2 == 0){
    jitter_amt <- c((-x/2):(x/2))[c((-x/2):(x/2))!=0]
  }
  dat_full_round_jittered[unique(jitter_list)[[i]], 2] <- dat_full_round_jittered[unique(jitter_list)[[i]], 2] + jitter_amt / ((abs(jitter_amt[1] - jitter_amt[2]))*10)
}

dat_full[,1] <- as.numeric(dat_full[,1])


#plot results - version 1, lines are thinner and longer
pdf(file = "yourfilename.pdf")
  plot(dat_full_round_jittered[,1:2], col = as.vector(dat_full_round_jittered[,3]), pch="_", cex=2, yaxt="n")
  par(yaxp= c(.5, 21.5, 21))
  axis(2)
  abline(c(.5,0))
  for (i in 1:21){
    abline(c(i+.5,0))
  }
  for (i in 1:10){
    abline(v= 250*(i-1) + 100 + as.numeric(max(genotypes_positioned[genotypes_positioned[,2]==i,3])))
  }
dev.off()

#plot results - version 2, lines are thicker and shorter
pdf(file = "yourfilename2.pdf")  
  plot(dat_full_round_jittered[,1:2], col = as.vector(dat_full_round_jittered[,3]), pch="-", cex=2.5, yaxt="n")
  par(yaxp= c(.5, 21.5, 21))
  axis(2)
  abline(c(.5,0))
  for (i in 1:21){
    abline(c(i+.5,0))
  }
  for (i in 1:10){
    abline(v= 250 *(i-1) + 100 + as.numeric(max(genotypes_positioned[genotypes_positioned[,2]==i,3])))
  }
dev.off()


#to edit these figures, I used a image editing software
#labeled the bottom axes chr 1-10
#labeled the y axes with elements
#cropped out extra lines on the side and blank portions on the top and bottom 



