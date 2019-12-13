
# This function obtains the data and adjust unbalance between 2 groups(optional)
# Genes with NA row or 0 frequency over 80% is deleted.
# Using variance of each gene, several top genes are selected.
# Returns data.frame of selected genes.


######################################################################
## Hyperparameter Setting #####################  topgene  ############
######################################################################
##------------------------------------------------------------------##
# Balancing data have 2 options
# Option 1 : Robust bootstrap from normal people
# Option 2 : Use some portion of cancer people only
##------------------------------------------------------------------##
top = 100           # How many genes will you eventually use?
balancing = FALSE   # Balancing number of cancer/normal 
balancing.opt = 1   # Option for balancing 
balance.ratio = 0.7 # Ratio of cancer people to total sample.
log = TRUE          # Print some messages.
plot = TRUE         # Print plots?
methodname = "RLE"  # Method of normalization / used in plot main names
######################################################################
######################################################################


# Row : gene / Column : People
topgene = function(genedata, top = 100, balancing = FALSE, balancing.opt = 1, balance.ratio = 0.7, log = TRUE, plot = TRUE, methodname = "KIHO"){
  
  #Printing Function Options
  cat("Selecting ", top, "genes", "from the dataset.\n", collapse="")
  cat("Option : ", collapse = "")
  if(!balancing){
    cat("No balancing.\n")
  }
  if(balancing){
    cat("Balancing with option", balancing.opt, "of ratio", balance.ratio, "\n\n", collapse="")
  }
  
  p = dim(genedata)[1]  # number of genes
  n = dim(genedata)[2]  # number of people

  #How many cancer / normal in the data
  name <- colnames(genedata)
  name_type <- substring(name, 14, 15)
    
  ##############################################
  #########     BALANCING PART    ##############
  ##############################################
  if(balancing){
    
    normal.index = (1:n)[name_type == "11"]   
    cancer.index = (1:n)[name_type == "01"]
    normal.n = length(normal.index)          #number of normal people
    cancer.n = length(cancer.index)          #number of cancer people
    
    if(balancing.opt == 1){
      boot.n = round((1/balance.ratio - 1) * cancer.n)  # how many sample should be bootstrapped from normal people
      bootindex = sample(normal.index, size = boot.n, replace = TRUE)
      bootstrapped = genedata[,bootindex]
      genedata = cbind(genedata[,cancer.index], bootstrapped)
    }
    
    if(balancing.opt == 2){
      select.n = round( (1/(1-balance.ratio)-1)*normal.n )     # how many sample should be selected from cancer people
      selectindex = sample(cancer.index, size = select.n, replace = FALSE)
      selecteddata = genedata[,selectindex]
      genedata = cbind(selecteddata, genedata[,normal.index])
    }
    
    # Just considering possible errors.
    if(!(balancing.opt %in% c(1,2))){
      cat("Error in balancing option number.\n")
      return()
    }
  }
  
  #Update number of samples 
  name <- colnames(genedata)
  name_type <- substring(name, 14, 15)
  n = dim(genedata)[2]  # number of people
  #table(name_type)
  
  # y : 0 if cancer 1 if normal
  y = ifelse(name_type =="01", 0, 1)
  #table(y)
  cat(paste("There are total ", n, "people with ", p , "genes in original data.\n", collapse=""))
  cat("Among them, ", table(name_type)["01"], "have cancer and ", table(name_type)["11"], "are normal.\n", collapse = "")
  
  
  #Selected index of gene (Keep updated through the code)
  selected = 1:p
  
  #Boxplot of gene expression  (by sample and by gene)
  if(plot){
    boxplot(genedata[,(1:20)*20], horizontal = TRUE, col = "red", xlab = "RLE normalized count", 
            main = paste("By sample", "(", methodname, ")", sep = "") )   #by sample
    boxplot(t(genedata[(1:20)*20,]), horizontal = TRUE, col = "red", xlab = "RLE normalized count", 
            main = paste("By gene", "(", methodname, ")", sep = ""))    #by gene
  }
  
  # Index of rows that contain NA (Maybe some error occurred during RLE)
  narow = (1:p)[is.na( apply(genedata, 1, sum) )]   # rows(genes) with NA value   
  data = genedata[setdiff(1:p, narow), ]  # select genes with no NA 
  selected = setdiff(1:p, narow)
  if(log){
    cat(paste("There are ", length(narow), "genes that have NA in the row.\n", collapse = ""))
  }
  
  # How many 0s for each gene
  zerofreq = apply(data, 1, function(x) sum(x==0))
  if(plot){
    plot(density(zerofreq), xlab = "# of people", main = paste("Frequency of 0 for each gene", "(", methodname, ")", sep = "") )
    abline(v = 0.8*n, col = "red", lty = 2) 
    ### Notice that the graph is bimodal, right subgraph starts approx. 0.8*n  (Red Dashed line)
  }
  
  #Filter genes with more than 80% of zeros
  manyzerogenes.n = sum(zerofreq >= 0.8*n)   #how many genes have genes less than 80% of 0 for total samples
  data2 = data[(zerofreq < 0.8*n), ]
  selected = selected[(zerofreq < 0.8*n)]
  if(log){
    cat(paste("There are ", manyzerogenes.n, "genes that have over 80% of 0 among people.\n", collapse = ""))
    cat(paste(length(selected),"genes are remaining.\n", collapse = ""))
  }

  
  #Caculate Variance and MAD(optional) for each genes
  genevar = apply(data2, 1, var)
  #genemad = apply(data2, 1, mad)
  if(plot){
    plot(density(genevar), xlab = "Variance of genes", main = paste("Variance of each gene(Density)", "(", methodname, ")", sep = "" ) )
    #plot(density(genemad), main = "MAD of each genes")
  }
  
  #TOP 100(hyperparameter) genes with big variance
  top.index = order(genevar, decreasing=TRUE)[1:top]
  data3 = data2[top.index, ]
  selected = selected[top.index]
  if(log){
    cat("Top ", top, "genes were selected successfully.\n", collapse = "" )
  }
  
  return(data3)
}
