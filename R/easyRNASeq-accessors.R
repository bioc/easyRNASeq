## If we have more than one readLength, return the one coresponding to the filename of the obj
## If we have one, return it
## In easyRNASeq, check that the names of the readLength are correct if lenhgth >1

##' Accessors for RNAseq class
##' 
##' 
##' These functions and generics define `accessors` (to get and set values) for
##' objects in the \pkg{easyRNASeq} package.
##' 
##' 
##' @aliases accessors chrSize chrSize,RNAseq-method
##' fileName<- fileName fileName,RNAseq-method fileName<-,RNAseq-method
##' geneModel geneModel<- geneModel,RNAseq-method geneModel<-,RNAseq-method
##' genomicAnnotation<- genomicAnnotation<-,RNAseq-method
##' genomicAnnotation genomicAnnotation,RNAseq-method
##' librarySize librarySize<- librarySize,RNAseq-method librarySize<-,RNAseq-method
##' organismName<- organismName organismName,RNAseq-method organismName<-,RNAseq-method 
##' readCounts<- readCounts readCounts,RNAseq-method readCounts<-,RNAseq-method
##' readCoverage<- readCoverage readCoverage,RNAseq-method readCoverage<-,RNAseq-method
##' readIslands<- readIslands readIslands,RNAseq-method readIslands<-,RNAseq-method
##' readLength<- readLength readLength,RNAseq-method readLength<-,RNAseq-method 
##' chrSize<- chrSize<-,RNAseq,integer-method chrSize<-,RNAseq,list-method 
##' seqnames,RNAseq-method
##' @name easyRNASeq accessors
##' @rdname easyRNASeq-accessors
##' @param obj An object derived from class \code{RNAseq}.
##' @param count The type of count you want to access,
##' 'genes','features','exons','transcripts' or 'islands'
##' @param summarization If count is set to genes, precise the type of
##' summarization, 'bestExons' or 'geneModels'
##' @param unique For the 'exons' count only. Should the counts returned be
##' unique for their identifier (i.e. the matrix row names)?
##' @param value The replacement value.
##' @usage genomicAnnotation(obj)
##' readCounts(obj,count=c("exons","features","genes","islands","transcripts"),
##' summarization=c("bestExons","geneModels"),unique=FALSE)
##' genomicAnnotation(obj) <- value
##' @return
##' Usually, the value of the corresponding slot, or other simple content
##' described on the help page of \code{easyRNASeq}.
##' @author Nicolas Delhomme
##' @keywords manip
##' @examples
##' 
##' rnaSeq<-new("RNAseq")
##' ##set organisme name of an RNAseq object
##' organismName(rnaSeq) <- "Dmelanogaster"
##' ##get organisme name of an RNAseq object
##' orgName<-organismName(rnaSeq)
##' 

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
              .Deprecated("datasource,AnnotParam-method",
                        msg="Getting the organism name is deprecated. Use an AnnotParam object instead and get its datasource.")
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
          definition=function(object){
            object@fileName
          })

setMethod(
          f="librarySize",
          signature="RNAseq",
          definition=function(obj){
            obj@librarySize
          })

##' @exportMethod seqnames
setMethod(
  f="seqnames",
  signature="RNAseq",
  definition=function(x){
    names(chrSize(x))
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
                 signature=c("RNAseq","list"),
                 definition=function(obj,value){

                   ## deprecated
                   .Deprecated("chrSize<-,RNAseq,numeric-methods",
                               msg="The use of the list for providing chromosome sizes has been deprecated. Use a named numeric vector instead.")
                                      
                   ## init
                   initialize(obj,chrSize=unlist(value))
                 })

setReplaceMethod(
                 f="chrSize",
                 signature=c("RNAseq","integer"),
                 definition=function(obj,value){
                   ## check
                   if(is.null(names(value))){
                     stop("We need a named vector or a named list for the chrSize slot!")
                   }
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
                     .Deprecated("datasource,AnnotParam-method",
                                 msg="Setting the organism name is deprecated. Use an AnnotParam object instead and set its datasource.")
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
