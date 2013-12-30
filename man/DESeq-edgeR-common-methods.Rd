\docType{methods}
\name{DESeq and edgeR common methods}
\alias{DESeq and edgeR common methods}
\alias{plotDispersionEstimates}
\alias{plotDispersionEstimates,CountDataSet-method}
\alias{plotDispersionEstimates,DGEList-method}
\title{DESeq and edgeR common methods}
\usage{
\S4method{plotDispersionEstimates}{CountDataSet}(obj, cond = character(1),
  log = "xy", ...)
}
\arguments{
  \item{obj}{An object of class
  \code{\linkS4class{CountDataSet}} or of class
  \code{\linkS4class{DGEList}}}

  \item{cond}{A character string describing the first
  condition; for
  \code{\linkS4class{CountDataSet}{CountDataSet}} obj
  only.}

  \item{log}{A character string passed onto
  \code{\link[graphics:plot.default]{plot.default}}; for
  \code{\linkS4class{CountDataSet}{CountDataSet}} obj
  only.}

  \item{...}{Additional plotting parameters; for
  \code{\linkS4class{CountDataSet}{CountDataSet}} obj
  only.}
}
\value{
none
}
\description{
\code{plotDispersionEstimates(obj,...)} extends the
\pkg{DESeq} and \pkg{edgeR} packages by offering the
functionality to plot the dispersion estimate as described
in their respective vignettes:
\code{\linkS4class{CountDataSet}{DESeq}} and
\code{\link[edgeR:edgeR-package]{edgeR}}.
}
\details{
\itemize{ \item \code{\linkS4class{CountDataSet}{DESeq}} A
character string describing the first condition, to be
provided as cond=value \item
\code{\link[edgeR:edgeR-package]{edgeR}} Unused, just for
compatibility.  }
}
\examples{
\dontrun{
	## edgeR
	## create the object
	dgeList <- DGEList(counts,group)
	## calculate the size factors
	dgeList <- calcNormFactors(dgeList)
	## plot them
	apply(combn(rownames(dgeList$samples),2),
		2,
		function(co,obj){plotNormalizationFactors(obj,co[1],co[2])},dgeList)
	## the dispersion estimates
	plotDispersionEstimates(obj)

	## DESeq
	## these are helper function for the DESeq package
	## refer to its vignette first
	cds <- newCountDataSet(countData,conditions)
	cds <- estimateSizeFactors(cds)
	cds <- estimateDispersions(cds)
	plotDispersionEstimates(cds,conditions[1])
	}
}
\author{
Nicolas Delhomme
}
\keyword{methods}

