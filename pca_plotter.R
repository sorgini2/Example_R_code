pca_plotter.f=function(numeric_matrix, char_matrix, char_matrix_col, dir, identifier,label_col) {
#This function will produce pca plots for the first 4 pc's and other useful pca stats from the following inputs
#numeric_matrix = The matrix of genotypes in numeric form (0, 1, 2 for the three genotypes).  A minimum of 10 markers is required.
#char_matrix = A table containing a list of unique identifiers for each sample in column 1, and remaining columns describing group identity.  Labels describing group identity should be no longer 18 characters.
#char_matrix_col = The column that contains the group identifiers by which you would like the samples grouped.  This will determine the colors and symbols used to group individuals.  There can be greater than 18 groupings.
#dir = name of the directory to which the data should be output.  This directory must exist already.
#identifier = string that will be used in the filename of the outputs
#label_col=the column that contains the labels for each point
  library(calibrate)
  group_names = sort(unique(as.character(char_matrix[,char_matrix_col])))
     
  colors = rainbow(6)  #Only 6 colours are truly distinguishable
  colors[2] = "#FFBF00FF" #And the yellow was too weak, so I got a darker one
  colors[6] = "black"
  colors = rep(colors, 3) #The maximum number of groupings is 18 - 6 colors x 3 symbols
  
  pch_list = c(rep(19,6), rep(17,6), rep(15,6)) #a list of point types for the plots that will ensure that when there are more than six individuals, colors get recycled and the point type changes (e.g. from circles to squares)
 pch_list[6]=1
  group_id = c() #Get a numeric vector of length = # of samples that indexes what group each sample belongs to
    for (i in 1:nrow(char_matrix)) {
      group_id[i] = which(group_names == char_matrix[i,char_matrix_col])
      }
  col_id = c()  #Get a numeric vector of length = # of samples that assigns a color to each sample according to group membership
    for (i in 1:nrow(char_matrix)) {
      col_id[i] = colors[group_id[i]]    
      }
  pch_id = c()
    for (i in 1:nrow(char_matrix)) {
      pch_id[i] = pch_list[group_id[i]]
      }
      
  #Scale the matrix
  pca_matrix = scale(numeric_matrix)
  #Set NAs to 0
  pca_matrix_noNA = apply(pca_matrix, 2, function (x) {ifelse(is.na(x), 0, x)})   

  #Perform PCA on pca_matrix_noNA
  eig = prcomp(pca_matrix_noNA, center = T)
  eigenvalues = (eig$sdev)^2
  scores = eig$x
  pc1_var = round((eigenvalues[1]/sum(eigenvalues)) * 100, 1)
  pc2_var = round((eigenvalues[2]/sum(eigenvalues)) * 100, 1)
  pc3_var = round((eigenvalues[3]/sum(eigenvalues)) * 100, 1)
  pc4_var = round((eigenvalues[4]/sum(eigenvalues)) * 100, 1)
  pc5_var = round((eigenvalues[5]/sum(eigenvalues)) * 100, 1)
  pc6_var = round((eigenvalues[6]/sum(eigenvalues)) * 100, 1)
  pc7_var = round((eigenvalues[7]/sum(eigenvalues)) * 100, 1)
  pc8_var = round((eigenvalues[8]/sum(eigenvalues)) * 100, 1)
  pc9_var = round((eigenvalues[9]/sum(eigenvalues)) * 100, 1)
  pc10_var = round((eigenvalues[10]/sum(eigenvalues)) * 100, 1)  
  #Plot the PC1 vs PC2
  png(paste(dir, identifier, "_", "pc1_vs_pc2.png", sep = ""), width = 10,units="in",res=200)
  layout(matrix(c(rep(1,18), rep(2,6)), 3, 8))
  par(mar = (c(5,5,4,0) + 0.1), pty = "m")
  plot(scores[,1], scores[,2], type = "n", xlab = paste("PC1 (", pc1_var, "%)", sep = ""), ylab = paste("PC2 (", pc2_var, "%)", sep = " "), cex.axis = 1.5, cex.lab = 1.5, main = "PC1 vs PC2", cex.main = 2)
  points(scores[,1], scores[,2], col = col_id, pch = pch_id, cex = 1.2)
  textxy(scores[,1], scores[,2],char_matrix[,label_col],cx = 0.5)
  par(mar = (c(5,0,4,0) + 0.1), pty = "m")
  plot(1:100, 1:100, type = "n", axes = F, xlab = "", ylab = "")
  rect(1, 0, 100, 100)
  points(rep(5, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), pch = pch_list[1:length(group_names)], col = colors[1:length(group_names)], cex = 1.1)
  text(rep(10, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), labels = gsub("_", " ", group_names), adj = 0)
  dev.off()
 
  
  #Plot the PC1 vs PC3
   png(paste(dir, identifier, "_", "pc1_vs_pc3.png", sep = ""), width = 10,units="in",res=200)
   layout(matrix(c(rep(1,18), rep(2,6)), 3, 8))
  par(mar = (c(5,5,4,0) + 0.1), pty = "m")
  plot(scores[,1], scores[,3], type = "n", xlab = paste("PC1 (", pc1_var, "%)", sep = ""), ylab = paste("PC3 (", pc3_var, "%)", sep = " "), cex.axis = 1.5, cex.lab = 1.5, main = "PC1 vs PC3", cex.main = 2)
  points(scores[,1], scores[,3], col = col_id, pch = pch_id, cex = 1.2)
  textxy(scores[,1], scores[,3],char_matrix[,label_col],cx=0.5)
  par(mar = (c(5,0,4,0) + 0.1), pty = "m")
  plot(1:100, 1:100, type = "n", axes = F, xlab = "", ylab = "")
  rect(1, 0, 100, 100)
  points(rep(5, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), pch = pch_list[1:length(group_names)], col = colors[1:length(group_names)], cex = 1.1)
  text(rep(10, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), labels = gsub("_", " ", group_names), adj = 0)
  dev.off()
 
  #Plot the PC1 vs PC4
    png(paste(dir, identifier, "_", "pc1_vs_pc4.png", sep = ""), width = 10,units="in",res=200)
  layout(matrix(c(rep(1,18), rep(2,6)), 3, 8))
  par(mar = (c(5,5,4,0) + 0.1), pty = "m")
  plot(scores[,1], scores[,4], type = "n", xlab = paste("PC1 (", pc1_var, "%)", sep = ""), ylab = paste("PC4 (", pc4_var, "%)", sep = " "), cex.axis = 1.5, cex.lab = 1.5, main = "PC1 vs PC4", cex.main = 2)
  points(scores[,1], scores[,4], col = col_id, pch = pch_id, cex = 1.2)
  textxy(scores[,1], scores[,4],char_matrix[,label_col],cx = 0.5)
  par(mar = (c(5,0,4,0) + 0.1), pty = "m")
  plot(1:100, 1:100, type = "n", axes = F, xlab = "", ylab = "")
  rect(1, 0, 100, 100)
  points(rep(5, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), pch = pch_list[1:length(group_names)], col = colors[1:length(group_names)], cex = 1.1)
  text(rep(10, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), labels = gsub("_", " ", group_names), adj = 0)
  dev.off()
 
   #Plot the PC1 vs PC5
  png(paste(dir, identifier, "_", "pc1_vs_pc5.png", sep = ""), width = 10,units="in",res=200)
    layout(matrix(c(rep(1,18), rep(2,6)), 3, 8))
  par(mar = (c(5,5,4,0) + 0.1), pty = "m")
  plot(scores[,1], scores[,5], type = "n", xlab = paste("PC1 (", pc1_var, "%)", sep = ""), ylab = paste("PC5 (", pc5_var, "%)", sep = " "), cex.axis = 1.5, cex.lab = 1.5, main = "PC1 vs PC5", cex.main = 2)
  points(scores[,1], scores[,5], col = col_id, pch = pch_id, cex = 1.2)
  textxy(scores[,1], scores[,5],char_matrix[,label_col],cx = 0.5)
  par(mar = (c(5,0,4,0) + 0.1), pty = "m")
  plot(1:100, 1:100, type = "n", axes = F, xlab = "", ylab = "")
  rect(1, 0, 100, 100)
  points(rep(5, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), pch = pch_list[1:length(group_names)], col = colors[1:length(group_names)], cex = 1.1)
  text(rep(10, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), labels = gsub("_", " ", group_names), adj = 0)
  dev.off()
 
    #Plot the PC1 vs PC6
   png(paste(dir, identifier, "_", "pc1_vs_pc6.png", sep = ""), width = 10,units="in",res=200)
  
  layout(matrix(c(rep(1,18), rep(2,6)), 3, 8))
  par(mar = (c(5,5,4,0) + 0.1), pty = "m")
  plot(scores[,1], scores[,6], type = "n", xlab = paste("PC1 (", pc1_var, "%)", sep = ""), ylab = paste("PC6 (", pc6_var, "%)", sep = " "), cex.axis = 1.5, cex.lab = 1.5, main = "PC1 vs PC6", cex.main = 2)
  points(scores[,1], scores[,6], col = col_id, pch = pch_id, cex = 1.2)
  textxy(scores[,1], scores[,6],char_matrix[,label_col],cx = 0.5)
  par(mar = (c(5,0,4,0) + 0.1), pty = "m")
  plot(1:100, 1:100, type = "n", axes = F, xlab = "", ylab = "")
  rect(1, 0, 100, 100)
  points(rep(5, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), pch = pch_list[1:length(group_names)], col = colors[1:length(group_names)], cex = 1.1)
  text(rep(10, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), labels = gsub("_", " ", group_names), adj = 0)
  dev.off()
 
     #Plot the PC1 vs PC7
   png(paste(dir, identifier, "_", "pc1_vs_pc7.png", sep = ""), width = 10,units="in",res=200)
   layout(matrix(c(rep(1,18), rep(2,6)), 3, 8))
  par(mar = (c(5,5,4,0) + 0.1), pty = "m")
  plot(scores[,1], scores[,7], type = "n", xlab = paste("PC1 (", pc1_var, "%)", sep = ""), ylab = paste("PC7 (", pc7_var, "%)", sep = " "), cex.axis = 1.5, cex.lab = 1.5, main = "PC1 vs PC7", cex.main = 2)
  points(scores[,1], scores[,7], col = col_id, pch = pch_id, cex = 1.2)
  textxy(scores[,1], scores[,7],char_matrix[,label_col],cx = 0.5)
  par(mar = (c(5,0,4,0) + 0.1), pty = "m")
  plot(1:100, 1:100, type = "n", axes = F, xlab = "", ylab = "")
  rect(1, 0, 100, 100)
  points(rep(5, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), pch = pch_list[1:length(group_names)], col = colors[1:length(group_names)], cex = 1.1)
  text(rep(10, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), labels = gsub("_", " ", group_names), adj = 0)
  dev.off()
 
 
     #Plot the PC1 vs PC8
   png(paste(dir, identifier, "_", "pc1_vs_pc8.png", sep = ""), width = 10,units="in",res=200)
   layout(matrix(c(rep(1,18), rep(2,6)), 3, 8))
  par(mar = (c(5,5,4,0) + 0.1), pty = "m")
  plot(scores[,1], scores[,8], type = "n", xlab = paste("PC1 (", pc1_var, "%)", sep = ""), ylab = paste("PC8 (", pc8_var, "%)", sep = " "), cex.axis = 1.5, cex.lab = 1.5, main = "PC1 vs PC8", cex.main = 2)
  points(scores[,1], scores[,8], col = col_id, pch = pch_id, cex = 1.2)
  textxy(scores[,1], scores[,8],char_matrix[,label_col],cx = 0.5)
  par(mar = (c(5,0,4,0) + 0.1), pty = "m")
  plot(1:100, 1:100, type = "n", axes = F, xlab = "", ylab = "")
  rect(1, 0, 100, 100)
  points(rep(5, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), pch = pch_list[1:length(group_names)], col = colors[1:length(group_names)], cex = 1.1)
  text(rep(10, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), labels = gsub("_", " ", group_names), adj = 0)
  dev.off()
 
     #Plot the PC1 vs PC9
  png(paste(dir, identifier, "_", "pc1_vs_pc9.png", sep = ""), width = 10,units="in",res=200)
    layout(matrix(c(rep(1,18), rep(2,6)), 3, 8))
  par(mar = (c(5,5,4,0) + 0.1), pty = "m")
  plot(scores[,1], scores[,9], type = "n", xlab = paste("PC1 (", pc1_var, "%)", sep = ""), ylab = paste("PC9 (", pc9_var, "%)", sep = " "), cex.axis = 1.5, cex.lab = 1.5, main = "PC1 vs PC9", cex.main = 2)
  points(scores[,1], scores[,9], col = col_id, pch = pch_id, cex = 1.2)
  textxy(scores[,1], scores[,9],char_matrix[,label_col],cx = 0.5)
  par(mar = (c(5,0,4,0) + 0.1), pty = "m")
  plot(1:100, 1:100, type = "n", axes = F, xlab = "", ylab = "")
  rect(1, 0, 100, 100)
  points(rep(5, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), pch = pch_list[1:length(group_names)], col = colors[1:length(group_names)], cex = 1.1)
  text(rep(10, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), labels = gsub("_", " ", group_names), adj = 0)
  dev.off()
 
     #Plot the PC1 vs PC10
  png(paste(dir, identifier, "_", "pc1_vs_pc10.png", sep = ""), width = 10,units="in",res=200)
    layout(matrix(c(rep(1,18), rep(2,6)), 3, 8))
  par(mar = (c(5,5,4,0) + 0.1), pty = "m")
  plot(scores[,1], scores[,10], type = "n", xlab = paste("PC1 (", pc1_var, "%)", sep = ""), ylab = paste("PC10 (", pc10_var, "%)", sep = " "), cex.axis = 1.5, cex.lab = 1.5, main = "PC1 vs PC10", cex.main = 2)
  points(scores[,1], scores[,10], col = col_id, pch = pch_id, cex = 1.2)
  textxy(scores[,1], scores[,10],char_matrix[,label_col],cx = 0.5)
  par(mar = (c(5,0,4,0) + 0.1), pty = "m")
  plot(1:100, 1:100, type = "n", axes = F, xlab = "", ylab = "")
  rect(1, 0, 100, 100)
  points(rep(5, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), pch = pch_list[1:length(group_names)], col = colors[1:length(group_names)], cex = 1.1)
  text(rep(10, length(group_names)), seq(95, 5, by = -(90/(length(group_names) - 1))), labels = gsub("_", " ", group_names), adj = 0)
  dev.off()
 
 
 
 
  #Boxplot for PC1 to PC10
  pdf(paste(dir, identifier, "_", "boxplots.pdf", sep = ""), width = 14)
  layout(matrix(1:10,2,5, byrow = T))
  for (i in 1:10) {
      scores_list = list()
      for (j in 1:length(group_names)) {
        scores_list[[j]] = scores[group_id == j, i]
        }
  par(mar = (c(3,8,2,1) + 0.1))
  boxplot(rev(scores_list), col = rev(colors[1:length(group_names)]), horizontal = T, yaxt = "n", main = paste("PC", i, sep = ""))
  axis(2, at = 1:length(group_names), labels = rev(gsub("_", " ", group_names)), las = 2, hadj = 1, cex.axis = 0.7)
  }
  dev.off()
}