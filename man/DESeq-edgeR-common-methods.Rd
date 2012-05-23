\name{DESeq and edgeR common methods}
\alias{DESeq and edgeR common methods}
\alias{plotDispersionEstimates}
\alias{plotDispersionEstimates,CountDataSet-method}
\alias{plotDispersionEstimates,DGEList-method}
\title{DESeq and edgeR common methods}
\arguments{
  \item{obj}{An object of class
  \code{\linkS4class{CountDataSet}} or of class
  \code{\linkS4class{DGEList}}}

  \item{\dots}{See details}
}
\value{
  none
}
\description{
  \code{plotDispersionEstimates(obj,...)} extends the
  \pkg{DESeq} and \pkg{edgeR} packages by offering the
  functionality to plot the dispersion estimate as
  described in their respective vignettes:
  \code{\linkS4class{CountDataSet}{DESeq}} and
  \code{\link[edgeR:edgeR-package]{edgeR}}.
}
\details{
  \itemize{ \item\code{\linkS4class{CountDataSet}{DESeq}} A
  character string describing the first condition, to be
  provided as cond=value
  \item\code{\link[edgeR:edgeR-package]{edgeR}} Unused,
  just for compatibility.  }
}
\examples{
\dontrun{
	## edgeR
	## create the object
	dgeList <- DGEList(counts,group)
	## calculate the sie factors
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

