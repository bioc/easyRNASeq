## internal functions

## conversion from aln, bam, GRanges to IRangesList
".extractIRangesList" <- function(obj,chr.sel=c()){
  ## the default is bam (easy, it has no class, just a list)
  aln.ranges <- switch(class(obj),
                       "AlignedRead" = {
                         if(any(is.na(position(obj)))){
                           stop("Your read file contains NA position, please filter for them. Check '?naPositionFilter'.")
                         }
                         ## no need to do it for the width, reads always have a positive width
                         split(IRanges(start=position(obj),width=width(obj)),as.character(chromosome(obj)))
                       },
                       "GappedAlignments" = split(ranges(obj),rname(obj)),
                       split(IRanges(start=obj$pos,width=obj$qwidth),obj$rname)
                       )
  if(!length(chr.sel)==0){
    aln.ranges <- aln.ranges[names(aln.ranges) %in% chr.sel,]
  }
  return(aln.ranges)
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
    val <- paste(",",paste(names(add.args[sel]),sapply(add.args,function(arg){ifelse(is.character(arg),paste("'",arg,"'",sep=""),arg)}),sep="=",collapse=", "))
  }
  return(val)
}

## to avoid reduce and strand errors (both exists in IRanges/intervals and GenomicRanges/genomeIntervals respectively)
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
                "RangedData"={unlist(values(genomicAnnotation(obj)[,sub("s$","",count)]))[,1]},
                stop(paste("No .getName functionality implemented for the class: ",class(genomicAnnotation(obj)),"!",sep=""))
                ))
}

