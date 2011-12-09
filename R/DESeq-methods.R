## this code was extracted from the DESeq vignette
## Simon Anders et al.

setMethod(
          f="fitInfo",
          signature="CountDataSet",
          definition=function(obj){
            obj@fitInfo
          })

setMethod(
          f="multivariateConditions",
          signature="CountDataSet",
          definition=function(obj){
            obj@multivariateConditions
          })

setMethod(
          f="plotDispersionEstimates",
          signature=c("CountDataSet"),
          definition=function(obj,cond=character(1),log="xy",...){
            
            ## check
            if(! class(obj) == "CountDataSet"){
              stop("This function: 'plotDispersionEstimates' only accepts a 'countDataSet' object.")
            }

            ## check the conditions
            if(multivariateConditions(obj)){
              ## TODO other steps are necessary before that to get the proper DESeq object at that point
              ## make sure we can read a data.frame as conditions, that the rownames are the files and that the estimates are calc with pooled-CR
              if(cond != "pooled"){
                stop("The provided condition can only have the value: 'pooled', as your conditions is multivariate")
              }
            } else {
              if(!cond %in% conditions(obj)){
                stop("The provided condition is not present in the 'conditions' slot of your object.")
              }
            }

            ## plot
            plot(
                 rowMeans( counts( obj, normalized=TRUE ) ),
                 fitInfo(obj)[[cond]]$perGeneDispEsts,
                 pch = '.', log=log,
                 xlab="gene mean normalized expression",
                 ylab="per gene dispersion estimate",
                 main = paste("Gene dispersion estimate for condition:",cond),
                 ...)
            xg <- 10^seq( -.5, 5, length.out=300 )
            lines( xg, fitInfo(obj)[[cond]]$dispFun( xg ), col="red" )
          })

