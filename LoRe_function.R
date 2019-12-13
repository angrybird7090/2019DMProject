
# This function conducts logistic regression on data, and test the model on testdata
# Output : ROC curve for predicted result on test data / Confusion Matrix with certain threshold (default 0.5)

########################################################################
## Hyperparameter Setting #########################  LoRe  #############
########################################################################
data = data.frame()     # Data with selected genes by topgene function
testdata = data.frame() # Data to test. Name of genes should be retained
threshold = 0.5         # Threshold value to print Confusion Matrix
plot = TRUE             # Plot ROC curve or not
methodname = "RLE"      # Method of normalization / used in plot main names
########################################################################
########################################################################


LoRe = function(data, testdata, threshold = 0.5, plot = TRUE, methodname = "KIHO"){
  library(pROC)
  library(caret)
  name <- colnames(data)
  name_type <- substring(name, 14, 15)
  # y : 0 if cancer 1 if normal
  y = ifelse(name_type =="01", 0, 1)
  
  top = dim(data)[1]
  data = as.data.frame(cbind(t(data),y))
  o <- train(as.factor(y) ~ ., data = data, method = "glm")
  
  # Predicting Test Data
  selected.genes = rownames(data)  # Which genes were selected by topgene function.
  test.y = ifelse(substring(colnames(testdata), 14, 15) =="01", 0, 1)
  test.y = as.factor(test.y)
  test.n = dim(testdata)[2]
  testdata = testdata[colnames(data)[1:top], ]
  
  testX = cbind(rep(1,test.n), t(testdata)) %*% o$finalModel$coefficients   #betax
  testX.prob = 1 - 1/(1+exp(testX))  # e^(betax) / (1 +  e^(betax))
  
  roc.test1 <- roc(test.y, as.vector(testX))
  if (log) {
    plot.roc(roc.test1, col="red", print.auc=TRUE, max.auc.polygon=TRUE, print.thres.pch=19, 
           print.thres.col = "red",  auc.polygon=TRUE, auc.polygon.col="#D1F2EB", 
           main = paste("LR ROC curve for testdata(", methodname,", Top", top, " genes)", sep ="" ))
    }   # 선 아래 면적에 대한 출력, 색상을 설정합니다. 
  roc.test.01result = ifelse(testX.prob > threshold, 1, 0)
  if (log){
    print(confusionMatrix(as.factor(roc.test.01result), as.factor(test.y))) }  #confusionMatrix
  return(confusionMatrix(as.factor(roc.test.01result), as.factor(test.y))$overall['Accuracy'])
}
