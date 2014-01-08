##' simpleRNASeq method
##' 
##' This function is a wrapper around the more low level functionalities of the
##' package. It is the simplest way to get a \code{\linkS4class{SummarizedExperiment}}
##' object from a set of bam files. \code{\linkS4class{SummarizedExperiment}} are
##' containers meant to hold any Next-Generation Sequencing experiment results and
##' metadata. The simpleRNASeq method replaces the 
##' \code{\link[easyRNASeq:easyRNASeq-easyRNASeq]{easyRNASeq}} function to 
##' simplify the usability. It does the following: \itemize{
##' \item use \code{\link[GenomicRanges:GRanges-class]{GenomicRanges}} 
##' for loading/pre-processing the data.
##' \item fetch \code{\link[easyRNASeq:fetchAnnotation]{annotations}}
##' depending on the selected parameters
##' \item extract the reads \code{\link[easyRNASeq:fetchCoverage]{coverage}} from
##' the provided file(s)
##' \item \code{\link[easyRNASeq:easyRNASeq-summarization-methods]{summarizes}} the
##' read counts according to the selected summarization
##' \item returns a \code{\linkS4class{SummarizedExperiment}} object.
##' }
##' 
##' @aliases simpleRNASeq simpleRNASeq,character-method
##' @rdname easyRNASeq-simpleRNASeq
##' @param bamFiles a \code{\linkS4class{BamFileList}} object
##' @param nnodes The number of CPU cores to use in parallel
##' @param param RnaSeqParam a \code{\linkS4class{RnaSeqParam}} object
##' that describres the RNA-Seq experimental setup.
##' @param verbose a logical to be report progress or not.
##' @return returns a \code{\linkS4class{SummarizedExperiment}} object.
##' @author Nicolas Delhomme
##' @seealso \code{\linkS4class{SummarizedExperiment}}
##' \code{\linkS4class{RnaSeqParam}}
##' @keywords methods
##' @examples
##' 
##'   \dontrun{
##'   TODO
##' }
##'
setMethod(f="simpleRNASeq",
          signature="character",
          definition=function(
          bamFiles=BamFileList(),
          param=RnaSeqParam(),
          nnodes=1,
          verbose=FALSE){
  
  ## create the output object
  sexp <- SummarizedExperiment(
    colData=DataFrame(
      FilePath=path(bamFiles),
      FileName=basename(names(bamFiles)),      
      row.names=basename(names(bamFiles))))
    
  ## validate the bams
  ## 1. the sequence info
  seqinfos <- lapply(bamFiles,seqinfo)
  if(length(seqinfos[[1]])==0){
    stop("It would seem some of your BAM files have no sequence information in their headers. Check these.")
  }
  exptData <- SimpleList(SeqInfo=seqinfos[[1]])
  if(!all(sapply(seqinfos,identical,exptData$SeqInfo))){
    stop("Your BAM files do not all have the same sequence informations. Check their headers.")
  }  
  exptData(sexp) <- exptData
  
  ## 2. paired or not and variable length or not  
  open(bamFiles)
  extracts <- lapply(lapply(bamFiles,scanBam,
                     param=ScanBamParam(what=c("qwidth","flag"),
                                        flag=scanBamFlag(isUnmappedQuery=FALSE))),
                     "[[",1)
  
  ## a. paired 
  colData(sexp)$Paired <- sapply(lapply(lapply(extracts,"[[","flag"),bamFlagTest,"isPaired"),any)
  
  ## b. width
  colData(sexp)$ReadLength <- sapply(lapply(lapply(extracts,"[[","qwidth"),range),
                                     function(rng){
                                       if(rng[1]==rng[2]){
                                         paste(rng[1],"bp",sep="")
                                         }else{
                                           paste(rng,"bp",sep="",collapse="-")
                                           }})
  close(bamFiles)
  
  ## 3. not a validation but add the read counts
  open(bamFiles)  
  colData(sexp)$TotalReads <- unlist(sapply(bamFiles,countBam)["records",])
  close(bamFiles)
  
  ## TODO determine strandedness from a number of non overlapping genes and check read orientation
  ## anything deviating from .5 would be assumed to be stranded
  
  ## TODO add the libSize
  ## LibSize=librarySize(obj),
  
  ## parallelize
  res <- parallelize(bamFiles,function(bamFile){

  ## open the bamFile
  open(bamFile)
    
  ## directly get the coverage?  
    
  ## stream the bam file
  .stream(bamFiles,yieldSize(param),verbose)
  
  ## if paired get the fragments
  
  ## get the coverage
  
  
  },nnodes)
})




