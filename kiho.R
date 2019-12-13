
# library----------------------------------------------------------------------
library(DESeq2)
library(pheatmap)
library(stats)
library(ggplot2)
library(knitr)
library(corrplot)
library(EDASeq)
library(rpart)
library(rpart.plot)
library(caret)

# read_data----------------------------------------------------------------------
read.data <- read.table("countdata.tsv", sep="\t", header = TRUE, row.names = 1)
read.data <- read.data[1:60483,] # 필요 없는 부분 제거
count.data <- round(2^read.data-1,0)

# row : gene
gene <- rownames(count.data)
tail <- substring(gene, 17, nchar(gene)) 
table(tail) # 염색체 확인

# col : name
name <- colnames(count.data)
table(substring(name, 14, 16)) # 사람들 type

# 사람들 type 01, 06, 11
name_type <- substring(name, 14, 15)
table(name_type)

cancer_index <- which(name_type=="01")
normal_index <- which(name_type=="11")

count.cancer <- count.data[, cancer_index]
count.normal <- count.data[, normal_index]


# fpkm
fpkm.read.data <- read.table("fpkmdata.tsv", sep="\t", header = TRUE, row.names = 1)
fpkm.read.data <- fpkm.read.data[,name]
fpkm.read.data <- fpkm.read.data[gene,]

fpkm.data <- 2^fpkm.read.data-1
fpkm.cancer <- fpkm.data[, cancer_index]
fpkm.normal <- fpkm.data[, normal_index]

# train set / test set-------------------------------------------------------------------------
nc <- length(cancer_index)
nn <- length(normal_index)

set.seed(2019)
train_cancer_index <- sample(1:nc, 731)
train_normal_index <- sample(1:nn, 75)

count.train <- cbind(count.cancer[,train_cancer_index], count.normal[,train_normal_index])
train_y <- c(rep(1,731), rep(0,75))
count.test <- cbind(count.cancer[,-train_cancer_index], count.normal[,-train_normal_index])
test_y <- c(rep(1,nc-731),rep(0,nn-75))

fpkm.train <- cbind(fpkm.cancer[,train_cancer_index], fpkm.normal[,train_normal_index])
fpkm.test <- cbind(fpkm.cancer[,-train_cancer_index], fpkm.normal[,-train_normal_index])

colData.train <- cbind(c(name[cancer_index][train_cancer_index], name[normal_index][train_normal_index]), train_y)
colData.train <- data.frame(colData.train)
colnames(colData.train) <- c("name", "group")
colData.test <- cbind(c(name[cancer_index][-train_cancer_index], name[normal_index][-train_normal_index]), test_y)
colData.test <- data.frame(colData.test)
colnames(colData.test) <- c("name", "group")

# normalization--------------------------------------------------------------------------
# tpm
tpm.train <- apply(fpkm.train, 2, function(x) x / sum(as.numeric(x)) * 10^6)
tpm.train.log <- log2(tpm.train + 1)

tpm.test <- apply(fpkm.test, 2, function(x) x / sum(as.numeric(x)) * 10^6)
tpm.test.log <- log2(tpm.test + 1)

colnames(tpm.train) <- colData.train$group
V <- apply(tpm.train, 1, var)
selectedGenes <- names(V[order(V, decreasing = T)][1:100])
pheatmap(tpm.train[selectedGenes,], scale = 'row', show_rownames = FALSE)

temp <- tpm.train[,c(1:10, 781:790)]
V <- apply(temp, 1, var)
selectedGenes <- names(V[order(V, decreasing = T)][1:1000])
pheatmap(temp[selectedGenes,], scale = 'row', show_rownames = FALSE)

#RLE
estimsf <- function (cts){
  # Compute the geometric mean
  geomMean <- function(x) exp(sum(log(x[x!=0]))/length(x[x!=0]))
  
  # Compute the geometric mean over the line
  gm.mean  <-  apply(cts, 1, geomMean)
  
  # Zero values are set to NA (avoid subsequentcdsdivision by 0)
  gm.mean[gm.mean == 0] <- NA
  
  # Divide each line by its corresponding geometric mean
  # sweep(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...)
  # MARGIN: 1 or 2 (line or columns)
  # STATS: a vector of length nrow(x) or ncol(x), depending on MARGIN
  # FUN: the function to be applied
  cts <- sweep(cts, 1, gm.mean, FUN="/")
  
  # Compute the median over the columns
  med <- apply(cts, 2, function(x) median(x[x!=0]))
  
  # Return the scaling factor
  return(med)
}
scaling <- estimsf(t(count.train))
rle.train = apply(count.train, 2, function(x) x / scaling)
rle.train.log = log2(rle.train + 1)
rle.train.log11 = log2(rle.train + 0.1)


scaling.test <- estimsf(t(count.test))
rle.test = apply(count.test, 2, function(x) x / scaling.test)
rle.test.log = log2(rle.test + 1)


#boxplot
par(mfrow=c(2,4))
boxplot(rle.train.log [,(1:20)*20],  col=train_y+1, cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="log rle")
boxplot(t(rle.train.log[(1:20)*20,]),  col=train_y+1, cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="log rle")

# log count
boxplot(read.data[,(1:20)*20],  col=train_y+1, cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="log count")
boxplot(t(read.data[(1:20)*20,]),  col=train_y+1, cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="log count")

boxplot(tpm.train.log[,(1:20)*20],  col=train_y+1, cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="log tpm")
boxplot(t(tpm.train.log[(1:20)*20,]),  col=train_y+1, cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="log tpm")

boxplot(fpkm.read.data[,(1:20)*20],  col=train_y+1, cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="log fpkm")
boxplot(t(fpkm.read.data[(1:20)*20,]),  col=train_y+1, cex.axis=0.7, 
        las=1, xlab="log2(counts)", horizontal=TRUE, main="log fpkm")

par(mfrow=c(1,1))

# pca --------------------------------------------------------------------------
M <- t(temp[selectedGenes,])
# transform the counts to log2 scale 
M <- log2(M + 1)
#compute PCA 
pcaResults <- prcomp(M)
coltemp <- colData.train[c(1:10, 781:790),]
plot(pcaResults)
summary(pcaResults)
screeplot(pcaResults, type="lines")

# correlation plot ------------------------------------------------------
correlationMatrix <- cor(temp)
kable(correlationMatrix,booktabs = TRUE)
corrplot(correlationMatrix, order = 'hclust')
corrplot(correlationMatrix, order = 'hclust', addrect = 2)
corrplot(correlationMatrix, order = 'hclust', addrect = 2, addCoef.col = 'white')
pheatmap(correlationMatrix)
pheatmap(correlationMatrix, annotation_col = coltemp)
pheatmap(correlationMatrix,  annotation_col = coltemp, cutree_cols = 2)

# Differential expression analysis--------------------------------------------
countData <- as.matrix(count.train)
designFormula <- "~ group"
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData.train, 
                              design = as.formula(designFormula))

print(dds)
rownames(dds)
colnames(dds)
counts(dds)
colData(dds)

dds.least.one <- dds[ rowSums(counts(dds)) > 1, ]
dds.deseq <- DESeq(dds.least.one)


#compute the contrast for the 'group' variable where 'CTRL' samples are used as the control group. 
DEresults = results(dds.deseq, contrast = c("group", '0', '1'))
#sort results by increasing p-value
DEresults <- DEresults[order(DEresults$pvalue),]
print(DEresults)
print(dds.deseq)

DESeq2::plotMA(object = dds.deseq, ylim = c(-7, 7))      
ggplot(data = as.data.frame(DEresults), aes(x = pvalue)) + geom_histogram(bins = 100)

countsNormalized <- counts(dds.deseq, normalized = TRUE)
selectedGenes <- names(sort(apply(countsNormalized, 1, var), decreasing = TRUE)[1:500])
plotPCA(countsNormalized[selectedGenes,], col = as.numeric(colData.train$group))

rld <- rlog(dds.deseq)
plotPCA(rld, ntop = 500, intgroup = 'group')


par(mfrow = c(1, 2))
plotRLE(countData, outline=FALSE, ylim=c(-4, 4), col=as.numeric(colData.train$group), main = 'Raw Counts')
plotRLE(counts(dds.deseq, normalized = TRUE), outline=FALSE, ylim=c(-4, 4), col = as.numeric(colData.train$group), ylim = c(-4,4), main = 'Normalized Counts')



# Functional Enrichment Analysis------------------------------------------------------------------
pheatmap(log2(M + 1), 
         annotation_col = colData[c(1:15, 797:806),], 
         show_rownames = FALSE, 
         scale = 'row', 
         cutree_cols = 2, 
         cutree_rows = 2)


# 망--------------------------------------------------------------------------------------------

#write.table(norm_train.log[which(tail==5),], "gene5.txt")

# 모두 0인 애들은 빼자
# va <- apply(traing_data,1,var)
# va[which(va==0)]
# me <- apply(traing_data,1,mean)
# max(me[which(va==0)])
# 
# train_data <- traing_data[-which(va==0),]
# new_gene <- gene[-which(va==0)]
  
# # 01 breast, 11 normal
# n1 <- sum(name_type=="01")
# n2 <- sum(name_type=="11")
# 
# # breast 데이터, normal 데이터
# breast <- un_dat_num[,which(name_type=="01")]
# normal <- un_dat_num[,which(name_type=="11")]
# 
# #breast mean과 normal mean
# breast_mean <- apply(breast, 1, mean)
# normal_mean <- apply(normal, 1, mean)
# 
# # breast, normal에서 var
# breast_var <- apply(breast, 1, var)
# normal_var <- apply(normal, 1, var)
# summary(breast_mean[which(breast_var==0 & normal_var==0)])
# summary(normal_mean[which(breast_var==0 & normal_var==0)])
# 
# table(apply(un_dat_num[,which(name_type=="06")], 1, mean)[which(breast_var==0 & normal_var==0)])
# a <- unname(data[which(apply(un_dat_num[,which(name_type=="06")], 1, mean)!=0 & breast_var==0 & normal_var==0),])
# a <- a[,-1]
# a[,which(name_type=="06")]
# 
# 
# index <- 1:nrow(un_dat_num)
# remove_index <- index[which(breast_var==0 & normal_var==0)]
# 
# # new data 만들고 그대로 다시
# new_data <- un_dat_num[-remove_index,]
# breast <- new_data[,which(name_type=="01")] 
# normal <- new_data[,which(name_type=="11")]
# 
# breast_mean <- apply(breast, 1, mean)
# normal_mean <- apply(normal, 1, mean)
# 
# breast_var <- apply(breast, 1, var)
# normal_var <- apply(normal, 1, var)
# 
# # t 통계량 밑에 들어갈 se
# se_bre_norm <- sqrt(breast_var/n1 + normal_var/n2)
# 
# # 차이만 가지고
# diff_bre_norm <- breast_mean - normal_mean
# summary(diff_bre_norm)
# hist(diff_bre_norm, breaks=1000)
# 
# # t 통계량 가지고
# t_bre_norm <- diff_bre_norm/se_bre_norm
# hist(t_bre_norm, breaks=100, probability = TRUE)
# lines(density(t_bre_norm), col=2, lwd=2)
# 
# 
# # quantile 들
# quantile(diff_bre_norm, 0.025)
# quantile(diff_bre_norm, 0.975)
# 
# quantile(t_bre_norm, 0.025)
# quantile(t_bre_norm, 0.975)
# 

# p-value 계산
pvalue <- rep(0,nrow(train_data))

for(i in 1:nrow(train_data)){
  r <- glm(train_y ~ train_data[i,], family = binomial)
  pvalue[i] <- summary(r)$coefficients[2,4]
}

hist(pvalue, breaks=500)
p_order <- order(pvalue)

# sub data 2000개 뽑기
sub_pvalue <- pvalue[p_order][1:2000]
sub_data <- train_data[p_order,][1:2000,]
sub_gene <- new_gene[p_order][1:2000]
hist(sub_pvalue, breaks=100)

# 염색체 확인
sub_tail <- substring(sub_gene, 17, nchar(sub_gene))
table(sub_tail)

# plot 그려보기
plot(sub_data[4,], sub_data[6,], col=-train_y+3, pch=train_y+2, lwd=-train_y+2)

# 다시 mean 들의 차이
sub_breast <- sub_data[,which(train_y==1)] 
sub_normal <- sub_data[,which(train_y==0)]
n1 <- ncol(sub_breast)
n2 <- ncol(sub_normal)

sub_breast_mean <- apply(sub_breast, 1, mean)
sub_normal_mean <- apply(sub_normal, 1, mean)

sub_breast_var <- apply(sub_breast, 1, var)
sub_normal_var <- apply(sub_normal, 1, var)

# t 통계량 밑에 들어갈 se
sub_se_bre_norm <- sqrt(sub_breast_var/n1 + sub_normal_var/n2)

# 차이만 가지고
sub_diff_bre_norm <- sub_breast_mean - sub_normal_mean
summary(sub_diff_bre_norm)
op <- par(mfrow=c(1,1))
hist(sub_diff_bre_norm, breaks=50)

# t 통계량 가지고
sub_t_bre_norm <- sub_diff_bre_norm/sub_se_bre_norm
hist(sub_t_bre_norm, probability = TRUE, breaks=50, xlab="t", main="Histogram of t statistic")
lines(density(sub_t_bre_norm), col=2, lwd=2)


# 뽑힌 것 중 바깥쪽
which(sub_t_bre_norm > 25)
which(sub_t_bre_norm < -30)

op <- par(mfrow=c(1,3))
plot(sub_data[703,], sub_data[886,], col=-train_y+3, pch=train_y+2, lwd=-train_y+2, xlab = sub_gene[703], ylab=sub_gene[886], main="Classification using genes with highest p-values")
legend("topright", legend=c("normal", "cancer"), pch = c(2,3), col = c(3,2))
plot(sub_data[905,], sub_data[870,], col=-train_y+3, pch=train_y+2, lwd=-train_y+2, xlab = sub_gene[905], ylab=sub_gene[870], main="Classification using genes with highest p-values")
legend("topright", legend=c("normal", "cancer"), pch = c(2,3), col = c(3,2))
plot(sub_data[905,], sub_data[886,], col=-train_y+3, pch=train_y+2, lwd=-train_y+2, xlab = sub_gene[905], ylab=sub_gene[886], main="Classification using genes with highest p-values")
legend("topright", legend=c("normal", "cancer"), pch = c(2,3), col = c(3,2))


# 뽑힌 것 중 안쪽
which(sub_t_bre_norm < 13 & sub_t_bre_norm > 11)
which(sub_t_bre_norm > -9 & sub_t_bre_norm < -6)

op <- par(mfrow=c(1,3))
plot(sub_data[921,], sub_data[998,], col=-train_y+3, pch=train_y+2, lwd=-train_y+2, xlab = sub_gene[921], ylab=sub_gene[998], main="Classification using genes with highest p-values")
legend("topright", legend=c("normal", "cancer"), pch = c(2,3), col = c(3,2))
plot(sub_data[999,], sub_data[772,], col=-train_y+3, pch=train_y+2, lwd=-train_y+2, xlab = sub_gene[999], ylab=sub_gene[772], main="Classification using genes with highest p-values")
legend("topright", legend=c("normal", "cancer"), pch = c(2,3), col = c(3,2))
plot(sub_data[772,], sub_data[921,], col=-train_y+3, pch=train_y+2, lwd=-train_y+2, xlab = sub_gene[772], ylab=sub_gene[921], main="Classification using genes with highest p-values")
legend("topright", legend=c("normal", "cancer"), pch = c(2,3), col = c(3,2))


# # vif 구하기
# stepwise_data <- step(glm(train_y ~ t(sub_data), family = binomial), steps=100, direction="backward")
# summary(stepwise_data)
# lambda_seq <-  10^seq(2, -2, by = -.1)
# cv.lambda <- cv.glmnet(t(sub_data), train_y, family="binomial", alpha=1, lambda = lambda_seq)
# summary(cv.lambda)
# a <- glm(train_y ~ t(sub_data[1:50,]), family = binomial)
# summary(a)



# svm
a <- svm(t(sub_data[1:1000,]), train_y,  scale = TRUE, type = "C-classification")

sub_test_data <- test_data[-which(va==0),][p_order,][1:1000,]
b <- predict(a, newdata = t(sub_test_data))
length(b)
hist(as.numeric(b)-1)
c <-as.numeric(b)-1
hist(c - test_y)

sum(abs(c-test_y))


# 망 2-----------------------------------------------------------------------------------------


dds <- DESeqDataSetFromMatrix(countData=train_data,
                              colData=colData,
                              design=~Group)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$Group <- relevel(dds$Group, ref=0)
dds <- DESeq(dds)
res <- results(dds, pAdjustMethod="BH")


eset <- exprs(as.matrix(data))
pset <- phenoData(as.matrix(data))





# 19.11.30

x <- train_data
colnames(x) <- train_y
conds = factor(train_y)
cds = newCountDataSet(x, conds) #컨디션 정보를 포함한 데이터 셋으로 저장
cds = estimateSizeFactors(x) #샘플간에 총 read 개수를 보정
cds = estimateDispersions(cds) #샘플간에 분산 및 오류를 예측
res = nbinomTest(cds, 'wt', 'mut') #샘플간 차등 발현 측정
resSig = res[res$padj<0.05] #보정된 p-value 가 0.05 인 것만 추출
write.csv(res, file='all_data.csv') #전체 데이터 파일에 쓰기
write.csv(resSig, file='sig_data.csv') 


#----------------------------------------------------------------------------


# tree---------------------------------------------------------

n <- 5000
a <- cbind(t(as.matrix(tpm.train[1:n,])),train_y)
dim(a)
colnames(a)[n+1] <- "y"
a <- as.data.frame(a)

default.ct = rpart(y~ ., data= a, method="class")
default.ct
prp(default.ct, type=1, extra=1, under = TRUE, split.font=1, varlen = -10,
    box.col = ifelse(default.ct$frame$var =="<leaf>", 'gray', 'white'))

default.ct.point.pred.train = predict(default.ct, a, type = "class")
confusionMatrix(default.ct.point.pred.train, factor(a$y))

full_dt=rpart(y~.,data=a, cp=0.1^20,method = "class")
printcp(full_dt)
plotcp(full_dt)
dt_prune <-prune(full_dt,cp=full_dt$cptable[which.min(full_dt$cptable[,"xerror"]),"CP"])
prp(dt_prune, type=1 , extra=1, under = T, split.font =2, varlen = -10, box.col = ifelse(dt_prune$frame$var == "<leaf>","gray","white"))

dt_prune.point.pred.train = predict(dt_prune, a, type="class")
confusionMatrix(dt_prune.point.pred.train, factor(a$y))

#---------------------------------------------------------



