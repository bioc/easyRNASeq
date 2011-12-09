## This code was extracted from the edgeR vignette
## Robinson et al.
setMethod(
          f="plotNormalizationFactors",
          signature=c("DGEList","character","character"),
          definition=function(obj=DGEList(),cond1=character(1),cond2=character(1)){
            
            ## check that the cond are in the samples
            if(nrow(obj$samples)==0){
              stop("Your 'DGEList' does not contain sample information!")
            }

            if(is.null(rownames(obj$samples))){
              stop("Your 'DGEList' samples data.frame has no row names!!")
            }

            if(!(all(c(cond1,cond2) %in% rownames(obj$samples)))){
              stop("The provided conditions do not exist in your 'DGEList' samples data.frame. Check the row names.")
            }
            
            ## raw data
            maPlot(obj$counts[, cond1], obj$counts[, cond2],
                   normalize = TRUE, pch = 19,
                   cex = 0.4, ylim = c(-8, 8),
                   main=paste("Smear plot",
                     rownames(obj$samples[cond1,]),
                     paste("(",obj$samples[cond1,"group"],")",sep=""),
                     "vs.\n",
                     rownames(obj$samples[cond2,]),
                     paste("(",obj$samples[cond2,"group"],")",sep=""),
                     "before normalization"))
            grid(col = "blue")
            abline(h = log2(obj$samples$norm.factors[match(cond1,rownames(obj$samples))]/obj$samples$norm.factors[match(cond2,rownames(obj$samples))]),
                   col = "red", lwd = 4)
            eff.libsize <- obj$samples$lib.size * obj$samples$norm.factors
            names(eff.libsize) <- rownames(obj$samples)
            ## norm data
            maPlot(obj$counts[, cond1]/eff.libsize[cond1], obj$counts[, cond2]/eff.libsize[cond2],
                   normalize = FALSE, pch = 19, cex = 0.4, ylim = c(-8, 8),
                   main=paste("Smear plot",
                     rownames(obj$samples[cond1,]),
                     paste("(",obj$samples[cond1,"group"],")",sep=""),
                     "vs.\n",
                     rownames(obj$samples[cond2,]),
                     paste("(",obj$samples[cond2,"group"],")",sep=""),
                     "after normalization"))
            grid(col = "blue") 
          })

setMethod(
          f="plotDispersionEstimates",
          signature=c("DGEList"),
          definition=function(obj=DGEList()){

            ## check that the cond are in the samples
            if(nrow(obj$samples)==0){
              stop("Your 'DGEList' does not contain sample information!")
            }

            ## check that the tagwise dispersion and the common dispersion were calculated
            if(is.null(obj$conc$conc.common)){
              stop("You need to estimate the common dispersion before using this function.")
            }

            if(is.null(obj$tagwise.dispersion)){
              stop("You need to estimate the tagwise dispersion before using this function")
            }

            ## plot
            plot(log2(obj$conc$conc.common), obj$tagwise.dispersion,
                 panel.first = grid(),ylab="tag dispersion",
                 xlab="common dispersion (log2)")
            abline(h = obj$common.dispersion, col = "dodgerblue", lwd = 3)
            legend("topright",col="dodgerblue",lwd=3,lty=1,bty="n","common dispersion")
          })

