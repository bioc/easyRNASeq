## this code was extracted from the DESeq vignette
## Simon Anders et al.

## not needed anymore, implemented in DESeq
## setMethod(
##           f="fitInfo",
##           signature="CountDataSet",
##           definition=function(obj){
##             obj@fitInfo
##           })


##' Extension for the DESeq package
##' 
##' \itemize{
##' \item{\code{multivariateConditions} is simply an accessor for the
##' \code{multivariateConditions} slot of a \code{\linkS4class{CountDataSet}}
##' object}
##' \item{\code{plotDispLSD} is a function silimar to 
##' \code{\link[DESeq:plotDispEsts]{plotDispEsts}}
##' that adds a density estimate as a colored heatmap from grey (few) to yellow
##' (many).}
##' }
##' 
##' 
##' @aliases multivariateConditions multivariateConditions,CountDataSet-method
##' plotDispLSD plotDispLSD,CountDataSet-method
##' @name DESeq additional methods
##' @rdname DESeq-methods
##' @param obj An object of class \code{\linkS4class{CountDataSet}}
##' @return \itemize{ \item{multivariateConditions returns a boolean describing
##' whether the data to analyze is multivariate or not }}
##' @author Nicolas Delhomme, Bastian Schiffthaler
##' @seealso \code{\linkS4class{CountDataSet}}
##' \code{\link[DESeq:plotDispEsts]{plotDispEsts}}
##' @keywords methods
##' @examples
##' 
##' 	\dontrun{
##' 	## these are helper function for the DESeq package
##' 	## refer to its vignette first
##' 	cds <- newCountDataSet(countData,conditions)
##' 	cds <- estimateSizeFactors(cds)
##' 	cds <- estimateDispersions(cds)
##' 	mVar <- multivariateConditions(cds)
##' 	}
##' 
setMethod(
          f="multivariateConditions",
          signature="CountDataSet",
          definition=function(obj){
            obj@multivariateConditions
          })

setMethod(
  f="plotDispLSD",
  signature="CountDataSet",
  definition=function(cds, name = NULL, ymin, 
                      linecol = "#00000080", xlab = "mean of normalized counts", 
                      ylab = "dispersion", log = "xy", cex = 0.45, ...) {
    px = rowMeans(counts(cds, normalized = TRUE))
    sel = (px > 0)
    px = px[sel]
    py = fitInfo(cds, name = name)$perGeneDispEsts[sel]
    if (missing(ymin)) 
      ymin = 10^floor(log10(min(py[py > 0], na.rm = TRUE)) - 
                        0.1)
    heatscatter(log10(px), log10(pmax(py, ymin)), xlab = xlab, ylab = ylab, 
                pch = ifelse(py < ymin, 6, 16), cexplot = cex, 
                xaxt='n', yaxt='n', ...)
    
    # Fix logged axes labels
    atx <- axTicks(1)
    aty <- axTicks(2)
    xlabels <- sapply(atx, function (i)
      as.expression(bquote(10^ .(i)))
    )
    ylabels <- sapply(aty, function (i)
      as.expression(bquote(10^ .(i)))
    )
    axis(1, at=atx, labels=xlabels)
    axis(2, at=aty, labels=ylabels)
    xg = 10^seq(-0.5, 5, length.out = 100)
    lines(log10(xg), log10(fitInfo(cds, name = name)$dispFun(xg)), col = linecol, 
          lwd = 4, lty = 1
    )
  })

##' DESeq and edgeR common methods
##' 
##' \code{plotDispersionEstimates(obj,...)} extends the \pkg{DESeq} and
##' \pkg{edgeR} packages by offering the functionality to plot the dispersion
##' estimate as described in their respective vignettes:
##' \code{\linkS4class{CountDataSet}{DESeq}} and
##' \code{\link[edgeR:edgeR-package]{edgeR}}.
##' 
##' \itemize{ \item \code{\linkS4class{CountDataSet}{DESeq}} A character string
##' describing the first condition, to be provided as cond=value
##' \item \code{\link[edgeR:edgeR-package]{edgeR}} Unused, just for
##' compatibility.  }
##' 
##' @aliases plotDispersionEstimates
##' plotDispersionEstimates,CountDataSet-method
##' plotDispersionEstimates,DGEList-method
##' @name DESeq and edgeR common methods 
##' @rdname DESeq-edgeR-common-methods 
##' @param obj An object of class \code{\linkS4class{CountDataSet}} or of class
##' \code{\linkS4class{DGEList}}
##' @param \dots See details
##' @return none
##' @author Nicolas Delhomme
##' @keywords methods
##' @examples
##' 
##' 	\dontrun{
##' 	## edgeR
##' 	## create the object
##' 	dgeList <- DGEList(counts,group)
##' 	## calculate the sie factors
##' 	dgeList <- calcNormFactors(dgeList)
##' 	## plot them
##' 	apply(combn(rownames(dgeList$samples),2),
##' 		2,
##' 		function(co,obj){plotNormalizationFactors(obj,co[1],co[2])},dgeList)
##' 	## the dispersion estimates
##' 	plotDispersionEstimates(obj)
##' 	
##' 	## DESeq
##' 	## these are helper function for the DESeq package
##' 	## refer to its vignette first
##' 	cds <- newCountDataSet(countData,conditions)
##' 	cds <- estimateSizeFactors(cds)
##' 	cds <- estimateDispersions(cds)
##' 	plotDispersionEstimates(cds,conditions[1])
##' 	}
##' 
setMethod(
          f="plotDispersionEstimates",
          signature=c("CountDataSet"),
          definition=function(obj,cond=character(1),log="xy",...){
            
            ## check
            if(! class(obj) == "CountDataSet"){
              stop("This function: 'plotDispersionEstimates' only accepts a 'countDataSet' object.")
            }

            ## check the conditions
            ## check if we are pooled, blind or per-condition
            stopifnot(length(cond)==1)

            ## should work without now
            ## if(multivariateConditions(obj)){
            ##   ## TODO other steps are necessary before that to get the proper DESeq object at that point
            ##   ## make sure we can read a data.frame as conditions, that the rownames are the files and that the estimates are calc with pooled-CR
            ##   if(cond != "pooled"){
            ##     stop("The provided condition can only have the value: 'pooled', as your conditions is multivariate")
            ##   }
            ## } else {
            ##   if(!cond %in% sub("disp_",names(fData(obj)))){
            ##     stop("The provided condition is not present in the 'conditions' slot of your object.")
            ##   }
            ## }

            ## plot
            plot(
                 rowMeans( counts( obj, normalized=TRUE ) ),
                 fitInfo(obj,cond)$perGeneDispEsts,
                 pch = '.', log=log,
                 xlab="gene mean normalized expression",
                 ylab="per gene dispersion estimate",
                 main = paste("Gene dispersion estimate for condition:",cond),
                 ...)
            xg <- 10^seq( -.5, 5, length.out=300 )
            lines( xg, fitInfo(obj,cond)$dispFun( xg ), col="red" )
          })

