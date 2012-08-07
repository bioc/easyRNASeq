\name{easyRNASeq annotation internal methods}
\alias{.geneModelAnnotation}
\alias{.getBmRange}
\alias{.getGffRange}
\alias{.getGtfRange}
\alias{.readGffGtf}
\alias{easyRNASeq annotation internal methods}
\title{Internal easyRNASeq annotation methods}
\usage{
  .getBmRange(organism = character(1), ...)
}
\arguments{
  \item{annotation.type}{describes the kind of annotation
  to keep the information from in a gtf or gff file. If set
  to NULL all the annotations are returned.}

  \item{fields}{added a parameter that allows defining the
  fields parsed from a gtf file. Still internal, but could
  easily be externalized.}

  \item{format}{describes the kind of annotation provided.
  One of gtf or gff.}

  \item{organism}{Organism name}

  \item{gAnnot}{a \code{\linkS4class{RangedData}} object
  containing exon annotations}

  \item{nbCore}{number of CPU cores to use}

  \item{filename}{filename that contains the annotations}

  \item{\dots}{Additional arguments, passed to more
  internal functions.}
}
\value{
  A \code{\linkS4class{RangedData}} containing the loaded
  or processed annotations.
}
\description{
  These are internal methods used to retrieve annotations
  tabularll \code{.getBmRange}Use
  \code{\link[biomaRt:useMart]{biomaRt}} to get exon
  annotations.  \code{.getGffRange}Use
  \code{\linkS4class{Genome_intervals_stranded}} to get
  annotation from a gff file.  \code{.getGtfRange}Use
  \code{\linkS4class{Genome_intervals_stranded}} to get
  annotation from a gtf file.
  \code{.geneModelAnnotation}Use the provided exon
  annotation to define gene models.  \code{.readGffGtf}Use
  \code{\linkS4class{Genome_intervals_stranded}} to get
  annotation from a gff or gtf file. It is called from
  \code{getGffRange} and \code{getGtfRange}.
}
\details{
  To use multicore machines more efficiently, the default
  parallel package will be used to parallelize the
  processing.
}
\author{
  Nicolas Delhomme
}
\keyword{internal}

