\name{DESeq methods}
\alias{fitInfo}
\alias{fitInfo,CountDataSet-method}
\alias{multivariateConditions}
\alias{multivariateConditions,CountDataSet-method}

\title{Extension for the DESeq package}
\description{
	\itemize{
	\item{\code{fitInfo}} is simply an accessor for the \code{fitInfo} slot of a \code{\linkS4class{CountDataSet}} object.
	\item{\code{multivariateConditions}} is simply an accessor for the \code{multivariateConditions} slot of a \code{\linkS4class{CountDataSet}} object.
	}
}

\usage{
fitInfo(obj)
multivariateConditions(obj)
}

\arguments{
  \item{obj}{An object of class \code{\linkS4class{CountDataSet}}}
  \item{\dots}{A character string describing the first condition, to be provided as cond=value}
}

\value{
	\itemize{
	\item{fitInfo} return the fit information
	\item{multivariateConditions} returns a boolean describing whether the data to analyze is multivariate or not
	}
}

\seealso{
\code{\linkS4class{CountDataSet}}
}

\examples{
	\dontrun{
	## these are helper function for the DESeq package
	## refer to its vignette first
	cds <- newCountDataSet(countData,conditions)
	cds <- estimateSizeFactors(cds)
	cds <- estimateDispersions(cds)
	fit.info <- fitInfo(cds)
	mVar <- multivariateConditions(cds)
	}
}

\author{Nicolas Delhomme}

\keyword{methods}

