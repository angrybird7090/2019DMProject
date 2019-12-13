
# This function conducts PCA on the data.
# Then, plots biplot of 2 PCs
# Finally, conduct LDA, QDA, NaiveBayes on PCAed data.

########################################################################
## Hyperparameter Setting #####################  pcaladaqda  ###########
########################################################################
axes = c(1,2)   # What component (PC) do you want to use for ggbiplot
                # axes should be vector of two integers
ldaqda = c(1,2,3) # 1 for lda, 2 for qda, 3 for NaiveBayes
                  # Multiple arguments such as c(1,2) is possible
precisionlda = 200 # Parameter for plotting LDA (more precise plotting)
methodname = "RLE"          # Method of normalization / used in plot main names
########################################################################
########################################################################

pcaldaqda = function(data, axes = c(1,2), ldaqda = c(1,2,3), precisionlda = 200, methodname = "KIHO"){
  library(ggbiplot)
  library(klaR)
  if(length(axes)!=2){
    cat("Need exactly 2 integers for axes.\n")
    return()
  }
  
  top = dim(data)[1]
  name <- colnames(data)
  name_type <- substring(name, 14, 15)
  # y : 0 if cancer 1 if normal
  y = ifelse(name_type =="01", 0, 1)
  
  pca = princomp(t(data))
  Label = factor(y, labels = c("Cancer", "Normal"))
  pcbiplot = ggbiplot(pca, choices = axes, var.axes=FALSE) + geom_point(aes(color = Label)) + 
    ggtitle(paste("Biplot of PC", axes[1], " and PC", axes[2], "(", methodname,", Top", top, " genes)", sep = ""))  #looks nice
  print(pcbiplot)
  print(screeplot(pca, type="lines", main = paste("Scree plot", "(", methodname,", Top", top, " genes)", sep = "") ))

  #Now LDA or QDA 
  s1 = pca$scores[,axes[1]]
  s2 = pca$scores[,axes[2]]
  partidata = as.data.frame(cbind(s2,s1))
  colnames(partidata) <- paste("PC", axes[c(2,1)], sep="")
  
  if(1 %in% ldaqda){
    partimat(as.factor(y) ~ ., data = partidata, method="lda",prec= precisionlda, 
             main = paste("Partition by LDA", "(", methodname,", Top", top, " genes)", sep = "") )
  }
  if(2 %in% ldaqda){
    partimat(as.factor(y) ~ ., data = partidata, method="qda",prec= precisionlda, 
             main = paste("Partition by QDA", "(", methodname,", Top", top, " genes)", sep = "")) 
  }
  if(3 %in% ldaqda){
    partimat(as.factor(y) ~ ., data = partidata, method="naiveBayes",prec= precisionlda, 
             main = paste("Partition by NaiveBayes", "(", methodname,", Top", top, " genes)", sep = "")) 
  }
}

