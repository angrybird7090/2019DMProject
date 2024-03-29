
# INPUTS

  # methodname: Method that was used to normalize the data (It is for main name of plots)
# topgene
  # genedata  : Normalized gene data
  # top       : How many genes will be selected eventually
  # balancing : Balancing number of samples
  # balancing.opt : Balancing option 1 or 2 ( 1:bootstrap(increase normal), 2:sampling(decrease cancer) )
  # balance.ratio : Preferring ratio of cancer/total
  # log       : Logical value for printing additional log messages or not
  # plot      : Logical value for printing plots or not
# pcaldaqda
  # data      : Data obtained by topgene function
  # axes      : What component (PC) do you want to use for ggbiplot  ex> c(1,2)
  #           : Should be vector of two integers
# ldaqda    : 1 for lda, 2 for qda, 3 for NaiveBayes
  #           : Multiple arguments such as c(1,2) is possible
  # precisionlda  : Parameter for plotting LDA (more precise plotting)
# LoRe
  # data      : Data obtained by topgene function
  # testdata  : Data to test. Name of genes should be retained
  # threshold : Threshold value to print Confusion Matrix



# OUTPUTS

# topgene   : 2 boxplot showing normalization / Frequency of 0 for each gene density graph / Variance of each gene(density)
# pcaldaqda : Scree plot / lda plot / qda plot/ naivebayes plot
# LoRe      : ROC curve of testdata / Confusion Matrix

# Adding Example
# Easy! Isn't it?

gene5  = read.table("gene5.txt")
testdata = read.table("test5.txt")

data11 = topgene(genedata = gene5, top = 100, balancing = TRUE, balancing.opt = 1, balance.ratio = 0.6, log = TRUE)
pcaldaqda(data = data11, testdata, axes = c(1,2), ldaqda = c(1,2,3), precisionlda = 200)
LoRe(data = data11, testdata = testdata, threshold = 0.5)   # This takes some time
myrf(data = data11, testdata = testdata, 0.5)

### With Cross Validation

estimated_test_error = cv(gene5, 5, top = 20, balancing.opt = 2, methodname = "Hyunsuk")

cv <- function(data, k, top = 100, balancing = TRUE, balancing.opt = 1, balance.ratio = 0.6, methodname = "KIHO"){
  
  n = ncol(data)
  name <- colnames(data)
  name_type <- substring(name, 14, 15)    
  normal.index = (1:n)[name_type == "11"]   
  cancer.index = (1:n)[name_type == "01"]
  normal.n = length(normal.index)          #number of normal people
  cancer.n = length(cancer.index)          #number of cancer people
  
  #cv_lab = sample(n, n, replace = F) %% k
  cv_lab_normal = sample(normal.n, normal.n, replace = F) %% k
  cv_lab_cancer = sample(cancer.n, cancer.n, replace = F) %% k
  
  # order : logistic regression, lda, qda, naive bayes, random forest
  cv_error = rep(0, 5)
  for (i_cv in 1:k){
    #w_val = which(cv_lab == (i_cv-1))
    w_val = c(normal.index[cv_lab_normal == (i_cv-1)], cancer.index[cv_lab_normal == (i_cv-1)])

    tr = data[,-w_val]
    val = data[,w_val]
    tr2 = topgene(genedata = tr, top = top, balancing = balancing, balancing.opt = balancing.opt, balance.ratio = balance.ratio, 
                  log = FALSE, plot = FALSE,  methodname = methodname)
    
    cv_error[1] = cv_error[1] + (1-LoRe(data = tr2, testdata = val, threshold = 0.5, log = FALSE, methodname = methodname))
    
    ldaqda_err_vec = pcaldaqda(data = tr2, val, axes = c(1,2), ldaqda = c(1,2,3), precisionlda = 200, plot.train = FALSE,  methodname = methodname)
    cv_error[2] = cv_error[2] + ldaqda_err_vec[1]
    cv_error[3] = cv_error[3] + ldaqda_err_vec[2]
    cv_error[4] = cv_error[4] + ldaqda_err_vec[3]
    
    cv_error[5] = cv_error[5] + (1-myrf(data = tr2, testdata = val, threshold = 0.5, log = FALSE,  methodname = methodname))
  }
  cv_error = cv_error/k
  
  print(cv_error)
}

#Also combining with MYRF, remove_topgene
library(dplyr)
data11 = topgene(genedata = gene5, top = 100, balancing = TRUE, balancing.opt = 2, balance.ratio = 0.6, log = FALSE, plot = FALSE) 
dat12 = data11 %>% remove_topgene(n = 20) 
LoRe(dat12, testdata)






#### ADDED BY HYUNSUK AFTER ADDING partimat2, pcaldaqda for test version, logistic regression using pca

#### ADDED FEATURES(INPUTS)
# All functions
#   - methodname  : Method that was used to normalize the data (It is for main name of plots)
#
# pcaldaqda
#   - plot.train  : Will you draw LDA/QDA plots for train data? (Default is FALSE, which implies drawing plots for test data)
#
# pcaLoRe
#   - pc.n        : Number of principle components you want to use for logistic regression
#                 : If given as vector, this function makes accuracy plot for each pc.n

gene5  = read.table("gene5.txt")
testdata = read.table("test5.txt")

data11 = topgene(genedata = gene5, top = 100, balancing = TRUE, balancing.opt = 1, balance.ratio = 0.6, log = TRUE, plot = TRUE, methodname = "KIHO")
# Train plot(as original pcaldaqda)
pcaldaqda(data = data11, testdata, axes = c(1,2), ldaqda = c(1,2,3), precisionlda = 200, plot.train = TRUE, methodname = "KIHO")
# Test plot(updated pcaldaqda)
pcaldaqda(data = data11, testdata, axes = c(1,2), ldaqda = c(1,2,3), precisionlda = 200, plot.train = FALSE, methodname = "KIHO")
LoRe(data = data11, testdata = testdata, threshold = 0.5)   # This takes some time
myrf(data = data11, testdata = testdata, 0.5)


pcaLoRe(data11, testdata, pc.n = 3, threshold = 0.5, plot = TRUE, methodname = "KIHO")
pcaLoRe(data11, testdata, pc.n = 1:10, threshold = 0.5, plot = TRUE, methodname = "KIHO")




