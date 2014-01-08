##' Extension for the edgeR package
##' 
##' This method extends the edgeR package by offering the functionality to plot
##' the effect of the normalization factor.
##' 
##' 
##' @aliases plotNormalizationFactors
##' plotNormalizationFactors,DGEList,character,character-method
##' @name edgeR additional methods
##' @rdname edgeR-methods
##' @param obj An object of class \code{\linkS4class{DGEList}}
##' @param cond1 A character string describing the first condition
##' @param cond2 A character string describing the second condition
##' @return none
##' @author Nicolas Delhomme
##' @keywords methods
##' @examples
##' 
##' 	\dontrun{
##' 	## create the object
##' 	dgeList <- DGEList(counts,group)
##' 	## calculate the sie factors
##' 	dgeList <- calcNormFactors(dgeList)
##' 	## plot them
##' 	apply(combn(rownames(dgeList$samples),2),
##' 		2,
##' 		function(co,obj){plotNormalizationFactors(obj,co[1],co[2])},dgeList)
##' 	}
##' 
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

            .Defunct(new="plotBCV",msg="This method is defunct due to change in
                     the edgeR DGEList object.")
          })

