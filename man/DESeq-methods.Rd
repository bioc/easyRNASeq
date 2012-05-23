\name{DESeq additional methods}
\alias{DESeq additional methods}
\alias{multivariateConditions}
\alias{multivariateConditions,CountDataSet-method}
\title{Extension for the DESeq package}
\arguments{
  \item{obj}{An object of class
  \code{\linkS4class{CountDataSet}}}
}
\value{
  \itemize{ \item{multivariateConditions returns a boolean
  describing whether the data to analyze is multivariate or
  not }}
}
\description{
  \itemize{ \item{\code{multivariateConditions} is simply
  an accessor for the \code{multivariateConditions} slot of
  a \code{\linkS4class{CountDataSet}} object.  }}
}
\examples{
\dontrun{
	## these are helper function for the DESeq package
	## refer to its vignette first
	cds <- newCountDataSet(countData,conditions)
	cds <- estimateSizeFactors(cds)
	cds <- estimateDispersions(cds)
	mVar <- multivariateConditions(cds)
	}
}
\author{
  Nicolas Delhomme
}
\seealso{
  \code{\linkS4class{CountDataSet}}
}
\keyword{methods}

