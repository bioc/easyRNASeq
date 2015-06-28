## this code was extracted from the DESeq vignette
## by Simon Anders et al.

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
##' \item{\code{plotDispersionEstimates} offers the functionality to plot
##' the dispersion estimate as described in the \pkg{DESeq} vignette.}
##' }
##'
##' @aliases multivariateConditions multivariateConditions,CountDataSet-method
##' plotDispLSD plotDispLSD,CountDataSet-method plotDispersionEstimates
##' plotDispersionEstimates,CountDataSet-method locfit lp newCountDataSet
##' @name DESeq additional methods
##' @rdname DESeq-methods
##' @param obj An object of class \code{\linkS4class{CountDataSet}}
##' @param cex standard \code{\link[graphics:plot.default]{plot.default}} parameter
##' @param cond A character string describing the first condition
##' @param linecol defines the line color
##' @param log A character string passed onto
##' \code{\link[graphics:plot.default]{plot.default}}
##' @param name argument passed to the \pkg{DESeq} \code{\link[DESeq:fitInfo]{fitInfo}}
##' function.
##' @param xlab standard \code{\link[graphics:plot.default]{plot.default}} parameter
##' @param ylab standard \code{\link[graphics:plot.default]{plot.default}} parameter
##' @param ymin numeric value defining the lower limit for the y axis
##' @param ... Additional plotting parameters
##' @usage multivariateConditions(obj)
##' plotDispLSD(obj, name = NULL, ymin,
##' linecol = "#00000080", xlab = "mean of normalized counts",
##' ylab = "dispersion", log = "xy", cex = 0.45, ...)
##' plotDispersionEstimates(obj,cond,log,...)
##' @return \itemize{
##' \item{\code{multivariateConditions} returns a boolean describing
##' whether the data to analyze is multivariate or not}
##' \item{\code{plotDispLSD} and \code{plotDispersionEstimates}} returns nothing
##' }
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
##' 	plotDispersionEstimates(cds,conditions[1])
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
  definition=function(obj, name = NULL, ymin,
                      linecol = "#00000080", xlab = "mean of normalized counts",
                      ylab = "dispersion", log = "xy", cex = 0.45, ...) {
    px = rowMeans(counts(obj, normalized = TRUE))
    sel = (px > 0)
    px = px[sel]
    py = fitInfo(obj, name = name)$perGeneDispEsts[sel]
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
    lines(log10(xg), log10(fitInfo(obj, name = name)$dispFun(xg)), col = linecol,
          lwd = 4, lty = 1
    )
  })

setMethod(
          f="plotDispersionEstimates",
          signature=c("CountDataSet"),
          definition=function(obj,cond=NULL,log="xy",...){

            ## check
            if(! class(obj) == "CountDataSet"){
              stop("This function: 'plotDispersionEstimates' only accepts a 'countDataSet' object.")
            }

            ## plot
            plot(
                 rowMeans( counts( obj, normalized=TRUE ) ),
                 fitInfo(obj,cond)$perGeneDispEsts,
                 pch = '.', log=log,
                 xlab="gene mean normalized expression",
                 ylab="per gene dispersion estimate",
                 main = paste("Gene dispersion estimate for",
                              ifelse(is.null(cond),"all conditions",
                                     paste("condition:",cond))),
                 ...)
            xg <- 10^seq( -.5, 5, length.out=300 )
            lines( xg, fitInfo(obj,cond)$dispFun( xg ), col="red" )
          })

