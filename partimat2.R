
# This function works almost same with partimat from klaR package
# It is upgrade version of it.
# This function draws partition plots for testdata


########################################################################
## Hyperparameter Setting #####################  partimat2  ############
########################################################################
# Others are same with original partimat function from klaR package
# testx         : Data to test (pc scores)
# testgrouping  : Real classification(answer of classification for test data)
########################################################################
########################################################################



partimat2 <- function(x, ...) 
  UseMethod("partimat2")


partimat2.default <- function(x, grouping, testx, testgrouping, stats = FALSE, method = "lda", prec = 100, 
                              nplots.vert, nplots.hor, main = "Partition Plot", name, mar,
                              plot.matrix = FALSE, plot.control = list(), ...){
  
  nvar <- ncol(x)
  if(nvar < 2) stop("at least 2 variables required")
  if(nlevels(grouping) < 2) stop("at least two classes required")
  nobs <- nrow(x)
  if(missing(name)) name <- colnames(x)
  # plot in scatterplot matrix
  
  ncomb <- round(0.5 * nvar * (nvar-1))
  if (missing(nplots.hor) && missing(nplots.vert)){
    nplots.hor<-ceiling(sqrt(ncomb))
    nplots.vert<-floor(sqrt(ncomb))
  }
  else if (missing(nplots.hor)) nplots.hor<-ceiling(ncomb/nplots.vert)
  else if (missing(nplots.vert)) nplots.vert<-ceiling(ncomb/nplots.hor)
  vars <- matrix(ncol=ncomb,nrow=2*nobs)
  varname <- matrix(ncol=ncomb,nrow=2)
  k <- 1
  for (i in 2:nvar)
    for (j in 1:(i-1))
    {
      vars[,k] <- c(x[,i], x[,j])
      varname[,k] <- c(name[i], name[j])
      k <- k + 1
    }
  
  if(missing(mar)) mar <- c(5.1, 4.1, 2.1, 1.1)
  opar <- par(mfrow = c(nplots.vert, nplots.hor), mar = mar, 
              oma = c(0, 0, !is.null(main), 0))
  on.exit(par(opar))
  
  sapply(1:ncomb, function(k) 
    drawparti2(grouping = grouping, testx = testx, testgrouping = testgrouping, stats = stats, x = vars[(1:nobs), k], 
               y = vars[(nobs+1):(2*nobs), k], method = method, 
               xlab = varname[1,k], ylab = varname[2,k], prec = prec, 
               legend.err = plot.matrix, plot.control = plot.control, ...)
  )
  par(mfrow=c(1,1))
  title(main = main, outer = TRUE)
  
  invisible()
}

partimat2.formula <- function(formula, data = NULL, ..., subset, na.action = na.fail) 
{
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval.parent(m$data))) 
    m$data <- as.data.frame(data)
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  grouping <- model.response(m)
  x <- model.matrix(Terms, m)
  xvars <- as.character(attr(Terms, "variables"))[-1]
  if ((yvar <- attr(Terms, "response")) > 0) 
    xvars <- xvars[-yvar]
  xlev <- if (length(xvars) > 0) {
    xlev <- lapply(m[xvars], levels)
    xlev[!sapply(xlev, is.null)]
  }
  xint <- match("(Intercept)", colnames(x), nomatch = 0)
  if (xint > 0) 
    x <- x[, -xint, drop = FALSE]
  res <- partimat2.default(x, grouping, ...)
  res$terms <- Terms
  cl <- match.call()
  cl[[1]] <- as.name("partimat")
  res$call <- cl
  res$contrasts <- attr(x, "contrasts")
  res$xlevels <- xlev
  attr(res, "na.message") <- attr(m, "na.message")
  if (!is.null(attr(m, "na.action"))) 
    res$na.action <- attr(m, "na.action")
  res
  invisible()
}

drawparti2 <- function(grouping, testx, testgrouping, x, y, stats = FALSE, method = "lda", prec = 100, 
                       xlab=NULL, ylab=NULL, col.correct = "black", col.wrong = "red", 
                       col.mean = "black", col.contour = "darkgrey", gs = as.character(grouping), testgs = as.character(testgrouping),
                       pch.mean = 19, cex.mean = 1.3, print.err = 0.7, legend.err = FALSE,
                       legend.bg = "white", imageplot = TRUE, image.colors = cm.colors(nc), 
                       plot.control = list(), ...){                       
  #grouping: class vec.
  #x: first data vec.
  #y: second data vec.
  #prec: nr. of hor/vert splits.
  
  success <- switch(method, 
                    rpart = requireNamespace("rpart") , 
                    naiveBayes = requireNamespace("e1071"))
  if(!is.null(success) && !success){
    message("For method 'rpart' the 'rpart' package is required, for method 'naiveBayes' the package 'e1071'.")
    return(NULL)
  }
  
  z <- switch(method,
              lda = lda(grouping ~ x + y,...),
              qda = qda(grouping ~ x + y,...),
              naiveBayes = e1071::naiveBayes(grouping~ x + y, 
                                             data = cbind.data.frame("grouping" = grouping, "x" = x, "y" = y), ...),
              stop("method not yet supported"))
  
  # Build a grid on the 2 coordinates
  xg <- seq(min(testx[,2]), max(testx[,2]), length = prec)
  yg <- seq(min(testx[,1]), max(testx[,1]), length = prec)
  grd <- expand.grid(x = xg, y = yg)
  # Calcultate posterior Probabilities on grid points
  temp <- switch(method,
                 lda = predict(z, grd,...)$post,
                 qda = predict(z, grd,...)$post,
                 naiveBayes = predict(z, grd , type="raw", ...),
                 stop("method not yet supported"))
  colnames(testx) <- c("y", "x")
  khead <- switch(method,
                  lda = predict(z, testx,...)$class,
                  qda = predict(z, testx,...)$class,
                  naiveBayes = predict(z, testx, ...),
                  stop("method not yet supported"))
  
  colorw <- testgrouping != khead
  err <- round(mean(colorw), 3)
  c00 = sum(testgrouping==0 & khead == 0)
  c01 = sum(testgrouping==0 & khead == 1)
  c10 = sum(testgrouping==1 & khead == 0)
  c11 = sum(testgrouping==1 & khead == 1)
  balance.err <- round(0.5*((c00/(c00+c01)) + (c11/(c10+c11))),3)

  #### If you want stats CURVE
  if(stats){
    library(caret)
    print(confusionMatrix(khead, testgrouping)$table)
    print(confusionMatrix(khead, testgrouping)$byClass)
  }
  
  color <- ifelse(colorw, col.wrong, col.correct)
  if(is.character(testgs) || is.factor(testgs)) testgs <- substr(testgs, 1, 1)
  
  nc <- ncol(temp)
  if(imageplot){
    do.call("image", c(list(xg, yg, matrix(apply(temp, 1, which.max), ncol = prec), 
                            main = NULL, col = image.colors, breaks = (0:nc) + .5, 
                            xlab = xlab, ylab = ylab), plot.control))
    do.call("points", c(list(testx[,2], testx[,1], pch = testgs, col = color), plot.control))
    box()
  }
  else 
    do.call("plot", c(list(testx[,2], testx[,1], pch = testgs, col = color, main = NULL, xlab = xlab, ylab = ylab), plot.control))
  #if((method=="lda") || (method=="qda")) 
   #points(z$means, pch = pch.mean, cex = cex.mean, col = col.mean)
  
  # For each class calculate the difference between prob. and max(prob) for other class,
  # so, the obs is assigned to class iff diff>0
  
  if(print.err){
    if(legend.err)
      legend(par("usr")[1], par("usr")[4], 
             legend = paste("Error:", err), bg = legend.bg, cex = print.err)
    else
      mtext(paste("app. balanced error rate:", 1-balance.err), 3, cex = print.err)    # use "app. error rate:", err instead of balance.err for Accuracy
  }
  
}