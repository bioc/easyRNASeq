#' Extension for the genomeIntervals package
#'
#' \describe{
#' \item{type}{ Another way to access the content of the gff type
#' column.  }}
#'
#' @aliases type,Genome_intervals-method
#' @name genomeIntervals additional methods
#' @rdname genomeIntervals-methods
#' @param x An object of class \code{\linkS4class{Genome_intervals}}
#' @usage \S4method{type}{Genome_intervals}(x)
#' @return \describe{
#' \item{type}{ The content of the
#' type column, usually a factor or a character vector } }
#' @examples
#'
#' library(curl)
#' library(genomeIntervals)
#' invisible(curl_download(paste0("https://github.com/UPSCb/UPSCb/raw/",
#'                                "master/tutorial/easyRNASeq/Dmel-mRNA-exon-r5.52.gff3.gz"),
#'                         "Dmel-mRNA-exon-r5.52.gff3.gz"))
#' annot<-readGff3("Dmel-mRNA-exon-r5.52.gff3.gz",quiet=TRUE)
#' type(annot)
#'
#' @author Nicolas Delhomme
#' @seealso
#' \itemize{
#' \item \code{\link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals
#' object}}
#' \item \code{\link[genomeIntervals:readGff3]{readGff3 function}}
#' }
# define a type accessor
setMethod(
          f="type",
          signature="Genome_intervals",
          definition=function(x){
            x$type
          })
