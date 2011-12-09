## TODO to allow different readLength, we can simply take a named list as readLength
## If we have more than one readLength, return the one coresponding to the filename of the obj
## If we have one, return it
## In easyRNASeq, check that the names of the readLength are correct if lenhgth >1

## getters
setMethod(
          f="genomicAnnotation",
          signature="RNAseq",
          definition=function(obj){
            obj@genomicAnnotation
          })

setMethod(
          f="readCoverage",
          signature="RNAseq",
          definition=function(obj){
            obj@readCoverage
          })

setMethod(
          f="readLength",
          signature="RNAseq",
          definition=function(obj){
            obj@readLength
          })

setMethod(
          f="chrSize",
          signature="RNAseq",
          definition=function(obj){
            obj@chrSize
          })

setMethod(
          f="readCounts",
          signature="RNAseq",
          definition=function(obj,
            count=c("exons","features","genes","islands","transcripts"),
            summarization=c("bestExons","geneModels"),
            unique=FALSE){
            
            ## If no count is given return the complete liste of count
            if(!length(count)==1){
              return(obj@readCounts)
            } else {
              .checkArguments("readCounts","count",count)
              return(switch(count,
                            "genes"= {                              
                              ## check that the summarization was provided
                              if(!length(summarization)==1){
                                stop(paste(
                                           "When the option 'count' is set to 'genes', one of the following summarization must be precised:",
                                           paste(eval(formals("readCounts")$summarization),collapse=",")
                                           )
                                     )
                              }
                              .checkArguments("readCounts","summarization",summarization)
                              ## return the matrix
                              return(switch(summarization,
                                            "bestExons"=obj@readCounts$genes$bestExons,
                                            "geneModels"=obj@readCounts$genes$geneModels
                                            ))
                            },
                            ## anything else returns the same:
                            {
                              tmp <- obj@readCounts[count][[1]]
                              if(unique){
                                tmp <- tmp[!duplicated(rownames(tmp)),,drop=FALSE]
                              }
                              return(tmp)
                            }
                            ## "islands"= obj@readCounts[count][[1]],                            
                            ## "transcripts"= obj@readCounts[count][[1]],
                            ))
            }
          })

setMethod(
          f="organismName",
          signature="RNAseq",
          definition=function(obj){
            obj@organismName
          })

setMethod(
          f="geneModel",
          signature="RNAseq",
          definition=function(obj){
            obj@geneModel
          })

setMethod(
          f="readIslands",
          signature="RNAseq",
          definition=function(obj){
            obj@readIslands
          })

setMethod(
          f="fileName",
          signature="RNAseq",
          definition=function(obj){
            obj@fileName
          })

setMethod(
          f="librarySize",
          signature="RNAseq",
          definition=function(obj){
            obj@librarySize
          })

## setters
setReplaceMethod(
                 f="genomicAnnotation",
                 signature="RNAseq",
                 definition=function(obj,value){
                   initialize(obj,genomicAnnotation=value)
                 })

setReplaceMethod(
                 f="readCoverage",
                 signature="RNAseq",
                 definition=function(obj,value){
                   initialize(obj,readCoverage=value)
                 })

setReplaceMethod(
                 f="readLength",
                 signature="RNAseq",
                 definition=function(obj,value){
                   if(is.numeric(value) & ! is.integer(value)){
                     value <- as.integer(value)
                     warning(paste("We expect an integer value and got a numeric one. Changing it to the integer: ", value, ".", sep=""))
                   }
                   initialize(obj,readLength=value)
                 })

setReplaceMethod(
                 f="chrSize",
                 signature="RNAseq",
                 definition=function(obj,value){
                   ## check
                   if(is.null(names(value))){
                     stop("We need a named vector or a named list for the chrSize slot!")
                   }
                   ## convert
                   if(is.vector(value)){
                     value <- as.list(value)
                   }
                   ## order
                   value[order(names(value))]
                   
                   ## init
                   initialize(obj,chrSize=value)
                 })

setReplaceMethod(
                 f="readCounts",
                 signature="RNAseq",
                 definition=function(obj,value){
                   initialize(obj,readCounts=value)
                 })

setReplaceMethod(
                 f="organismName",
                 signature="RNAseq",
                 definition=function(obj,value){
                   initialize(obj,organismName=value)
                 })

setReplaceMethod(
                 f="geneModel",
                 signature="RNAseq",
                 definition=function(obj,value){
                   initialize(obj,geneModel=value)
                 })

setReplaceMethod(
                 f="readIslands",
                 signature="RNAseq",
                 definition=function(obj,value){
                   initialize(obj,readIslands=value)
                 })

setReplaceMethod(
                 f="fileName",
                 signature="RNAseq",
                 definition=function(obj,value){
                   initialize(obj,fileName=value)                   
                 })

setReplaceMethod(
                 f="librarySize",
                 signature="RNAseq",
                 definition=function(obj,value){
                   initialize(obj,librarySize=value)                   
                 })
