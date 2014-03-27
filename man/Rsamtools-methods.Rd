\docType{methods}
\name{validate,BamFile-method}
\alias{validate}
\alias{validate,BamFile-method}
\alias{validate,BamFileList-method}
\title{Extension of the Rsamtools package}
\usage{
\S4method{validate}{BamFile}(obj, header = TRUE, cross.validation = TRUE)
}
\arguments{
  \item{obj}{An object of the \code{\linkS4class{BamFile}}
  or \code{\linkS4class{BamFileList}} class}

  \item{header}{a boolean to (de)activate the check for a
  BAM header}

  \item{cross.validation}{a boolean - only valid for
  \code{\linkS4class{BamFileList}} objects - to
  (de)activate the cross validation of all the BAM files
  header}
}
\value{
\code{validate} returns invisibly a vector of boolean.
Fails anyway if any file is missing.
}
\description{
Describes extensions to the Rsamtools package. \itemize{
\item For \code{\linkS4class{BamFile}} and
\code{\linkS4class{BamFileList}} objects: \itemize{ \item
\code{validate} validates a \code{\linkS4class{BamFile}} or
\code{\linkS4class{BamFileList}} object. } }
}
\details{
\code{validate} checks whether the BAM file exists and if a
BAI index is present.
}
\examples{
filenames <- dir(path=system.file("extdata",
      		                  package="RnaSeqTutorial"),
					       pattern="[A,C,T,G]{6}\\\\.bam$",
  				       full.names=TRUE)

bfl <-BamFileList(filenames,index=filenames)

validate(bfl)
}
\author{
Nicolas Delhomme
}
\seealso{
\itemize{ \item \code{\linkS4class{BamFile}} \item
\code{\linkS4class{BamFileList}} }
}
\keyword{methods}

