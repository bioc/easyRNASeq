% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/easyRNASeq-simpleRNASeq.R
\name{simpleRNASeq,BamFileList,RnaSeqParam-method}
\alias{simpleRNASeq,BamFileList,RnaSeqParam-method}
\alias{simpleRNASeq}
\title{simpleRNASeq method}
\usage{
\S4method{simpleRNASeq}{BamFileList,RnaSeqParam}(
  bamFiles = BamFileList(),
  param = RnaSeqParam(),
  nnodes = 1,
  verbose = TRUE,
  override = FALSE
)
}
\arguments{
\item{bamFiles}{a \code{\linkS4class{BamFileList}} object}

\item{param}{RnaSeqParam a \code{\linkS4class{RnaSeqParam}} object
that describes the RNA-Seq experimental setup.}

\item{nnodes}{The number of CPU cores to use in parallel}

\item{verbose}{a logical to be report progress or not.}

\item{override}{Should the provided parameters override the detected ones}
}
\value{
returns a \code{\linkS4class{RangedSummarizedExperiment}} object.
}
\description{
This function is a wrapper around the more low level functionalities of the
package. It is the simplest way to get a \code{\linkS4class{RangedSummarizedExperiment}}
object from a set of bam files. \code{\linkS4class{RangedSummarizedExperiment}} are
containers meant to hold any Next-Generation Sequencing experiment results and
metadata. The simpleRNASeq method replaces the
\code{\link[easyRNASeq:easyRNASeq-easyRNASeq]{easyRNASeq}} function to
simplify the usability. It does the following:
\itemize{
\item use \code{\link[GenomicAlignments:GAlignments-class]{GenomicAlignments}}
for reading/pre-processing the BAM files.
\item get the \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{annotations}}
depending on the selected parameters
\item calculate the coverage from the provided file(s)
\item \code{\link[easyRNASeq:easyRNASeq-summarization-methods]{summarizes}} the
read counts according to the selected summarization
\item returns a \code{\linkS4class{RangedSummarizedExperiment}} object.
}
}
\examples{

# the data
tdir <- tutorialData()
annot <- fetchData("Drosophila_melanogaster.BDGP5.77.with-chr.gtf.gz")

 # create the BamFileList, get the BAM and BAI index files from the Bioc cache
 filenames <- dir(tdir,pattern="[A,T].*\\\\.bam$",full.names=TRUE)
 indexnames <- sapply(paste0(sub(".*_","",basename(filenames)),".bai"),fetchData)
 bamFiles <- getBamFileList(filenames,indexnames)

  # create the AnnotParam
  annotParam <- AnnotParam(annot,type="gtf")

  # create the RnaSeqParam
  rnaSeqParam <- RnaSeqParam(annotParam=annotParam)

  # get a RangedSummarizedExperiment containing the counts table
  sexp <- simpleRNASeq(
    bamFiles=bamFiles,
    param=rnaSeqParam,
    verbose=TRUE
  )

  # get the counts
  assays(sexp)$exons

}
\seealso{
\itemize{
\item{For the input:
\itemize{
\item \code{\linkS4class{AnnotParam}}
\item \code{\linkS4class{BamParam}}
\item \code{\linkS4class{RnaSeqParam}}
}}
\item{For the output:
\code{\linkS4class{RangedSummarizedExperiment}}
}
\item{For related functions:
\itemize{
\item \code{\linkS4class{BamFile}}
\item \code{\linkS4class{BamFileList}}
\code{\link[easyRNASeq:easyRNASeq-BamFileList]{getBamFileList}}
}
}}
}
\author{
Nicolas Delhomme
}
\keyword{methods}
