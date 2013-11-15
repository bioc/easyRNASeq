\name{GenomicRanges additional methods}
\alias{GenomicRanges additional methods}
\alias{colnames}
\alias{colnames,GRangesList-method}
\alias{colnames,GenomicRanges-method}
\title{Extension of the GenomicRanges package}
\arguments{
  \item{obj}{An object of the \code{\linkS4class{GRanges}}
  or \code{\linkS4class{GRangesList}} class}
}
\value{
  A vector of column names.
}
\description{
  Return the column name of a \code{\linkS4class{GRanges}}
  or \code{\linkS4class{GRangesList}} object.
}
\details{
  It returns the actual column names of the elementMetadata
  slot of the \code{\linkS4class{GRanges}} or
  \code{\linkS4class{GRangesList}} object. The
  elementMetadata contains a \code{\linkS4class{DataFrame}}
  object used to store additional information provided by
  the user, such as exon ID in our case.
}
\examples{
\dontrun{
	## an example of a RangedData annotation
	gAnnot <- RangedData(
                     IRanges(
                             start=c(10,30,100),
                             end=c(21,53,123)),
                          space=c("chr01","chr01","chr02"),
                          strand=c("+","+","-"),
                          transcript=c("trA1","trA2","trB"),
                          gene=c("gA","gA","gB"),
                          exon=c("e1","e2","e3"),
                          universe = "Hs19"
                          )

	## an example of a GRangesList annotation
	grngs <- as(gAnnot,"GRanges")

	## accessing the colnames
	colnames(grngs)

	## creating a GRangesList
	grngsList<-split(grngs,seqnames(grngs))

	## accessing the colnames
	colnames(grngsList)
	}
}
\author{
  Nicolas Delhomme
}
\seealso{
  \code{\linkS4class{DataFrame}}
  \code{\linkS4class{GRanges}}
  \code{\linkS4class{GRangesList}}
}
\keyword{methods}

