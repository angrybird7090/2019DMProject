
# This function conducts random forest on data, and test the model on testdata

# Input : data      : Data obtained by topgene function
#         testdata  : Data to test. Name of genes should be retained
#         threshold : Threshold value to print Confusion Matrix
#         plot      : Plot ROC curve or not
#         methodname: Method of normalization / used in plot main names
# Output : ROC curve for predicted result on test data / Confusion Matrix with certain threshold (default 0.5)

myrf = function(data, testdata, threshold = 0.5, plot = TRUE, methodname = "KIHO"){
  library(pROC)
  library(caret)
  library(randomForest)

  name <- colnames(data)
  name_type <- substring(name, 14, 15)
  # y : 0 if cancer 1 if normal
  y = ifelse(name_type =="01", 0, 1)

  top = dim(data)[1]
  data = as.data.frame(cbind(t(data),y))
  rf<-randomForest(as.factor(y) ~ ., data=data,ntree=500)
  # Predicting Test Data
  selected.genes = rownames(data)  # Which genes were selected by topgene function.
  test.y = ifelse(substring(colnames(testdata), 14, 15) =="01", 0, 1)
  test.y = as.factor(test.y)
  test.n = dim(testdata)[2]
  testdata = testdata[colnames(data)[1:top], ]
  test.x = t(testdata) 
  result.predicted.prob <- predict(rf, test.x, type="vote")
  
  # Draw ROC curve.
  result.roc <- roc(test.y, result.predicted.prob[,2])
  if(plot){
    #plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
    plot.roc(result.roc, col="red", print.auc=TRUE, max.auc.polygon=TRUE, print.thres.pch=19, 
             print.thres.col = "red",  auc.polygon=TRUE, auc.polygon.col="#D1F2EB", 
             main = paste("RF ROC curve for testdata(", methodname,", Top", top, " genes)", sep ="" ))
  }
  pred.y = predict(rf, test.x)
  # Confusion Matrix
  roc.test.result = ifelse(result.predicted.prob[,1] > threshold, 0, 1)
  if (log) {print(confusionMatrix(as.factor(roc.test.result), as.factor(test.y)))}
  return(confusionMatrix(as.factor(roc.test.result), as.factor(test.y))$overall['Accuracy'])
}
