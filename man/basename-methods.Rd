% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basename-methods.R
\docType{methods}
\name{basename methods}
\alias{basename methods}
\alias{basename}
\alias{basename,BamFile-method}
\alias{basename,BamFileList-method}
\title{Extend the basename function to display Rsamtools BamFile class basename}
\usage{
\S4method{basename}{BamFile}(path)
}
\arguments{
\item{path}{an object of class \code{\linkS4class{BamFile}} or \code{\linkS4class{BamFileList}}}
}
\description{
Display the basename of the bam file represented by a \code{\linkS4class{BamFile}} object.
}
\section{Methods}{
 \describe{ \item{list("signature(object = \"BamFile\")")}{
Display the basename of the bam file linked to by a
\code{\linkS4class{BamFile}} object.  } }
}

\keyword{methods}
