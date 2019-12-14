
# Code to save plots easily


########################################################################
## Hyperparameter Setting ################  topgene.plots.save  ########
########################################################################

#Other parameters are parameters for topgene function
#Added paramters
data.array = c("data1", "data2")  # Array of data names (strings) to conduct topgene, 
                                  # in the order of plotname.mehtods(c("Count","RLE", "TPM", "FPKM"))
filename.prefix = "PREFIX"        # Your choice for prefix of filename of saving plots

########################################################################
########################################################################


### Example Code ###

topgene.plots.save(data.list = list(countdata, rledata, tpmdata, fpkmdata), filename.prefix = "WHATPREFIX")





### Function Code ###

normalize.methods = c("Not", "RLE", "TPM", "FPKM")
plotname.methods = c("Count","RLE", "TPM", "FPKM")
data.list = list()  # List of dataframes!!  ex> list(countdata, rledata, tpmdata, fpkmdata)


topgene.plots.save = function(data.list, filename.prefix = "PREFIX", top = 100, balancing = TRUE, balancing.opt = 1, balance.ratio = 0.5, log = FALSE, 
                              normalize.methods = c("Not", "RLE", "TPM", "FPKM"), plotname.methods = c("Count","RLE", "TPM", "FPKM")){
  
  for( i in 1:length(normalize.methods) ){
    
    cat("\nPlotting data \t\t:", data.array[i],"\n")
    cat("Normalized method \t:", plotname.methods[i],"\n")
    
    png(paste(filename.prefix, "_", plotname.methods[i], ".png", sep = ""), width = 800, height = 600)
    par(mfrow=c(2,2))
    data.count <- topgene(genedata = data.list[[i]], top = top, balancing = balancing, 
                          balancing.opt = balancing.opt, balance.ratio = balance.ratio, log = log, plot = TRUE, methodname = normalize.methods[i])
    cat("-----------------------------------------------------------------------\n")
    dev.off()
  }
  
  
}



