## TODO update the returned values in the Roxygen doc
##' Internal methods of RNAseq objects
##' 
##' These are generic internal methods:
##' \itemize{
##' \item .catn Just some pretty printing.
##' \item .checkArguments check that the provided argument match one
##' of the formal definition of the function. Stop if not.
##' \item .convertToUCSC
##' convert chromosome names to UCSC compliant ones.
##' \item .extractIRangesList extract an IRanges object from an AlignedRead or a
##' GAlignments object or a list returned by reading a bam file with
##' Rsamtools. It returns a list containing the IRangesList and library size.
##' \item .getArguments For a given function returns the arguments
##' passed as part of the \dots{} that match that function formals.
##' \item .getName Get the genomicAnnotation object names. Necessary to deal
##' with the different possible annotation object: \code{RangedData} or
##' \code{GRangesList}.
##' \item .list.files check the arguments passed through the \dots to select only the valid ones.
##' \item .normalizationDispatcher a function to dispatch
##' the normalization depending on the 'outputFormat' chosen by the user.
##' \item reduce Allow proper dispatch between the
##' \link[intervals:Intervals_virtual-class]{intervals} and the
##' \link[GenomicRanges:GRanges-class]{GenomicRanges} reduce function
##' \item strand Allow proper dispatch between the
##' \link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals} and
##' the \link[GenomicRanges:GRanges-class]{GenomicRanges} strand function
##' \item strand<- Allow proper dispatch between the
##' \link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals} and
##' the \link[GenomicRanges:GRanges-class]{GenomicRanges} strand replace
##' function
##' }
##' 
##' @aliases .catn .checkArguments .convertToUCSC 
##' .extractIRangesList .getArguments .getName
##' .list.files
##' .normalizationDispatcher reduce reduce,RNAseq-method strand
##' strand,RNAseq-method strand<- strand<-,RNAseq-method
##' @name easyRNASeq internal methods
##' @rdname easyRNASeq-internal-methods
##' @param arg The argument name to check for.
##' @param chr.names The chromosome names, as a character vector, to be
##' converted to UCSC ones
##' @param chr.sel A list of chromosome to restrict the IRanges spaces
##' returned.
##' @param fun The name of the function
##' @param obj An RNAseq object, or for the 'normalizationDispatcher',
##' depending on the type: a CountDataSet, a DGEList, a matrix, or an RNAseq
##' object respectively
##' @param organism The organism name
##' @param type character string specifying the type of object
##' (normalizationDispatcher)
##' @param value the appropriate strand object (strand and strand<-) or the
##' provided argument value (checkArguments)
##' @param x an object of the
##' \link[GenomicRanges:GRanges-class]{GenomicRanges},
##' \link[intervals:Intervals_virtual-class]{intervals} or
##' \link[genomeIntervals:Genome_intervals_stranded-class]{genomeIntervals}
##' package
##' @param ... For \code{.getArguments} a list of named parameters to be
##' matched against a function formal definition. For \code{.catn}, the values
##' to be printed. For \code{.list.files}, the additional parameters to be filtered
##' for the list.files function.
##' @return
##' \item{argString}{a character string representing these arguments
##' and their value that matched those defined in the formal definition of the
##' function}
##' \item{convertedChrNames}{a converted vector of chromosome names}
##' \item{i.range}{an IRange object} \item{names}{The
##' annotation names, i.e. a combination of exon, feature, transcript and gene}
##' \item{normalized.counts}{Depending on the type, a CountDataSet, a DGEList,
##' a NumericList, or NULL respectively}
##' @author Nicolas Delhomme
##' @keywords internal
## conversion from aln, bam, GRanges to IRangesList
## returns as well the library size
".extractIRangesList" <- function(obj,chr.sel=c()){
  ## the default is bam (easy, it has no class, just a list)
  ## and only the aligned reads are kept at this point
  
  ## pre-filter
  if(!length(chr.sel)==0){
    obj <- switch(class(obj),
                  "AlignedRead" = {
                    obj[chromosome(obj) %in% chr.sel,]
                  },
                  "GAlignments" = {
                    obj[seqnames(obj) %in% chr.sel,]
                  },
                  {
                    sel <- obj$rname %in% chr.sel
                    lapply(obj,"[",sel)
                  }
                  )      
  }

  ## TODO check if countBam is not faster for both bam and GAlignments
  
  ## get the ranges and lib sizes
  return(switch(class(obj),
                "AlignedRead" = {
                  if(any(is.na(position(obj)))){
                    stop("Your read file contains NA position, please filter for them. Check '?naPositionFilter'.")
                  }
                  
                  ## warn for multiple mapping
                  if(sum(duplicated(sread))>0){
                    warning("Your alignment file potentially contains multi-mapping reads. This would bias the counting.")
                  }
                  
                  ## no need to do it for the width, reads always have a positive width
                  list(rng.list=split(IRanges(start=position(obj),width=width(obj)),as.character(chromosome(obj))),
                       lib.size=length(obj))
                },
                "GAlignments" = {
                  ## we want only the mapped regions
                  grnglist <- grglist(obj,drop.D.range=TRUE)
                  
                  ## check if NH tag are present
                  ## we check only for the first...
                  if(! "NH" %in% colnames(mcols(grnglist[[1]]))){
                    warning("Your alignment file misses the NH tag. It may contains multi-mapping reads, which would bias the counting.")
                  } else {                  
                    if(any(as.data.frame(grnglist)$NH)>1){
                      warning("Your alignment file potentially contains multi-mapping reads. This would bias the counting.")
                    }
                  }
                  list(rng.list=split(unlist(ranges(grnglist),use.names=FALSE),unlist(seqnames(grnglist),use.names=FALSE)),
                       lib.size=length(unique(names(obj))))
                },
                {
                  if(any(obj$tag$NH)>1){
                    warning("Your alignment file potentially contains multi-mapping reads. This would bias the counting.")
                  }
                  list(rng.list=split(IRanges(start=obj$pos,width=obj$qwidth),obj$rname),
                     lib.size=length(obj$rname))
                }
                ))
}

## keep only the valid names
## TODO offer the possibility for the user to give a function to resolve that
## TODO check why chr.names is a list and not a c()
".convertToUCSC" <- function(chr.names=list(),
                            organism=c("Dmelanogaster","Hsapiens","Mmusculus","Rnorvegicus"),
                            custom.map=data.frame()){

  if(is.factor(chr.names)){
    chr.names <- as.character(chr.names)
  }
  
  if(organism!="custom"){
    ## TODO preprocess the names
    ## e.g. remove the .fa extension

    ## check that we do not already have valid names
    if(length(grep("chr",chr.names))==length(chr.names) & length(grep("MT|dmel_mitochondrion_genome",chr.names))==0){
      return(chr.names)			
    }
  }
  
  return(switch(tolower(organism),
                "dmelanogaster"={
                  paste(
                        "chr",
                        sub("dmel_mitochondrion_genome","M",chr.names),
                        sep="")
                },
                "hsapiens"={
                  paste(
                        "chr",
                        sub("MT","M",chr.names),
                        sep="")
                },
                "mmusculus"={
                  paste(
                        "chr",
                        sub("MT","M",chr.names),
                        sep="")
                },
                "rnorvegicus"={
                  paste(
                        "chr",
                        sub("MT","M",chr.names),
                        sep="")
                },
                "custom"={

                  ## sanity checks
                  if(!is.data.frame(custom.map)){
                    stop("Your custom map is not a data.frame as expected")                    
                  }

                  if(length(colnames(custom.map)) != 2){
                    stop("Your custom map does not have the expected number (2) of columns")
                  }
                  
                  if(!all(colnames(custom.map) %in% c("from","to"))){
                    stop("Your custom map does not follow the column names convention. They should be named 'from' and 'to'.")
                  }

                  ## convert to character
                  if(is.factor(custom.map$from) | is.factor(custom.map$to)){
                    custom.map <- data.frame(apply(custom.map,2,as.character),stringsAsFactors=FALSE)
                  }
                    
                  ## in case there are already valid
                  if(length(na.omit(match(chr.names,custom.map$to))) != length(custom.map$to)){
                    ## if not convert them
                    sel <- match(chr.names,custom.map$from)
                    if(all(is.na(sel))){
                      stop(paste("Your custom map does not match any chromosome names in the list:",paste(chr.names,collapse=", ")))
                    }
                    chr.names[!is.na(sel)] <- na.omit(custom.map[sel,"to"])
                    if(any(is.na(sel))){
                      warning(paste("Your custom map does not define a mapping for the following chromosome names:",paste(chr.names[is.na(sel)],collapse=", ")))
                    }
                  }
                  chr.names
                },
                {
                  warning(paste("No function implemented to convert the names for the organism: ",organism,".\n",
                                "Available ones so far exists for: ",paste(eval(formals(".convertToUCSC")[["organism"]]),collapse=", "),".\n",
                                "If you need to provide your own mapping, use the 'custom' argument.",sep=""))
                  chr.names
                }))
}

## check arguments
".checkArguments" <- function(fun,arg,value){
  args <- eval(formals(fun)[[arg]])

  ## first stop if too many values
  if(length(value)!=1){
    stop(paste(
               "Please select a value for the argument: '",arg,
               "' among the supported arguments: ",
               paste(args,collapse=", "),sep=""))
  }

  ## second stop if not valid
  if(!value %in% args){
    stop(paste(
               "The given ",arg," argument:",
               value,
               "is not part of the supported arguments:",
               paste(args,collapse=", ")))
  }
}

## get arguments
".getArguments" <- function(fun,...){
  ## the possible args
  args <- names(eval(formals(fun)))

  ## the provided args
  add.args <- list(...)

  ## filter
  sel <- names(add.args) %in% args

  ## report them as a string
  val<-""
  if(sum(sel)>0){
    val <- paste(",",paste(names(add.args[sel]),sapply(add.args[sel],function(arg){ifelse(is.character(arg),paste("'",arg,"'",sep=""),arg)}),sep="=",collapse=", "))
  }
  return(val)
}

## to avoid reduce and strand errors (both exists in IRanges/intervals and GenomicRanges/genomeIntervals respectively)
##' @exportMethod strand
setMethod(
          f="strand",
          signature="RNAseq",
          definition=function(x){
            if(extends(class(x),"Genome_intervals_stranded")){
              genomeIntervals::strand(x)
            } else {
              GenomicRanges::strand(x)
            }
          })
##' @exportMethod reduce
setMethod(
          f="reduce",
          signature="RNAseq",
          definition=function(x){
            if(extends(class(x),"Intervals_virtual")){
              genomeIntervals::reduce(x)
            } else {
              GenomicRanges::reduce(x)
            }
          })
##' @exportMethod strand<-
setReplaceMethod(
                 f="strand",
                 signature="RNAseq",
                 definition=function(x,value){
                   if(extends(class(x),"Genome_intervals_stranded")){
                     genomeIntervals::strand(x) <- value
                   } else {
                     GenomicRanges::strand(x)<-value
                   }
                   return(x)
                 })

## pretty printing
".catn" <- function(...){
  cat(...,"\n")
}

## get the names (exon ID, transcript ID, ...) from either a RangedData or GRangesList
## NB this is possible, because we order the read object (read coverage) according to the annotations!!!!
".getName" <- function(obj,count){
  return(switch(class(genomicAnnotation(obj)),
                "GRangesList"={unlist(
                                      lapply(
                                             lapply(
                                                    lapply(genomicAnnotation(obj),
                                                           values),
                                                    "[",match(sub("s$","",count),colnames(genomicAnnotation(obj)))),
                                             unlist),
                                      use.names=FALSE)},
                "RangedData"={
                  unlist(values(genomicAnnotation(obj)[,sub("s$","",count)])[,1],use.names=FALSE)
                  },
                stop(paste("No .getName functionality implemented for the class: ",class(genomicAnnotation(obj)),"!",sep=""))
                ))
}

## to allow list.files additional parameters
.list.files <- function(path,pattern,...){
  dots <- list(...)
  dots <- dots[names(dots) %in% names(formals(fun=list.files))]
  dots <- dots[!names(dots) %in% c("path","pattern","full.names")]
  l.args <- c(list(path=path,pattern=pattern,full.names=TRUE),dots)
  
  eval(parse(text=
             gsub("[\\]","\\\\\\\\",
                  paste("list.files(",
                        paste(names(l.args),paste("",l.args,"",sep="'"),
                              sep='=',collapse=","),")",sep=""))
             ))
}

