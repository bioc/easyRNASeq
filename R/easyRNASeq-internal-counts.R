##' Summarization of reads count per gene
##' 
##' Called by the method \code{geneCounts}. Having performed the
##' \code{exonCounts} is a pre-requisite.  \itemize{
##' \item.bestExonSummarization Identify the exon showing the highest coverage.
##' \item.geneModelSummarization Sum the coverage values of the synthetic exons
##' constituting a gene model.  }
##' 
##' 
##' @aliases .bestExonSummarization .geneModelSummarization
##' @name easyRNASeq summarization internal methods
##' @rdname easyRNASeq-summarization-internal-methods
##' @param obj An object derived from class \code{\linkS4class{RNAseq}}
##' @return A vector containing read counts
##' @author Nicolas Delhomme
##' @keywords internal
## TODO move the Rd from easyRNASeq-internal-methods.R to here for the count internal methods
".bestExonSummarization" <- function(obj){
	## check if we have the exon summary already
	if(is.null(readCounts(obj,'exons'))){
		exonCounts(obj)
	}
	eCounts<-data.frame(                            
                            .getName(obj,"genes"),
                            .getName(obj,"exons"),
                            readCounts(obj,'exons'))
	colnames(eCounts) = c("gene","exon","counts")
	eCounts$cov<-eCounts$counts/width(genomicAnnotation(obj))
	eCounts<-eCounts[order(eCounts$gene,eCounts$cov,decreasing=TRUE),]
	sel<-!duplicated(eCounts$gene)
	gCounts <- eCounts[sel,"counts"]
	names(gCounts) <- eCounts[sel,"exon"]
	return(gCounts)
}

".geneModelSummarization" <- function(obj){
  
  ## select for the proper spaces
  gm.sel <- na.omit(match(names(readCoverage(obj)),names(geneModel(obj))))
  
  gAgg <- stats:::aggregate(
                            as.integer(
                                       unlist(
                                              aggregate(
                                                        readCoverage(obj)[names(geneModel(obj))[gm.sel]],
                                                        ranges(geneModel(obj))[gm.sel],
                                                        sum)
                                              )
                                       ),
                            list(gene=geneModel(obj)[gm.sel]$gene),
                            sum)
  gCounts <- ceiling(gAgg[,2])
  names(gCounts)<-gAgg$gene
  return(gCounts)
}

".doBasicCount" <- function(obj){

  ##check genomicAnnotation
  if(length(genomicAnnotation(obj))==0){
    stop('The genomicAnnotation slot is empty')
  }
  ##check readCoverage
  if(length(readCoverage(obj))==0){
    stop('The readCoverage slot is empty')
  }
  
  ## not used by supposedly faster than aggregate
  ##RleViewsList(rleList=trackCoverage, rangesList=exons)
            
  ## counts
  counts=as.integer(unlist(aggregate(
    readCoverage(obj)[match(names(genomicAnnotation(obj)),names(readCoverage(obj)))],
    ranges(obj),
    sum
    )))

  ## return value
  return(ceiling(counts))
}

".extendCountList" <- function(
                              cList,values,
                              type=c("exons","features","transcripts","genes","islands"),
                              subType=character(1),
                              filename=character(1)){
	
  ## check the type
  types <- eval(formals(".extendCountList")$type)
  if(!type %in% types){
    stop(paste(
               "The given type:",
               type,
               "is not part of the supported types:",
               paste(types,collapse=", ")))
  }

  ## check the summarization
  if(length(subType)>1){
    warning(paste(
                  "We can only handle one 'subType'! Re-setting the 'subType' from:",
                  paste(subType,sep=", "),"to ''"))
    subType=character(1)
  }
  
  ## Special treatment for Islands counts
  if(type == 'islands'){
    if(length(cList)>0){
      filename=paste(filename,type,sep='.')
      eval(parse(text=paste("cList$",type,"$filename<- values",sep="")))
    }
    else{
      eval(parse(text=paste("cList$",type,"<- values",sep="")))
    }
    
  }
  ##check if the object is initialized or not
  else if(length(cList)>0){
    values <- as.matrix(values)
#    colnames(values) <- paste(filename,type,sep=".")
    ## Set clist with value for different conditions
    check <- function(list,add){
      if(length(list)>0){
        ## we need a selector
        sel <- colnames(add) %in% colnames(list)
        if(all(sel)){
          list <- add
        } else { 
          if(any(sel)){
            list[,colnames(add)[sel]] <- add[,sel]
          }
          
          list <- cbind(list,add[,!sel])
          rownames(list) <- c(rownames(list),names(add)[!sel])
        }
      }
      else{
        list <- add
      }
      return(list)
    }
    ## Use check with cList$type$subType and values
    eval(parse(text=paste("cList$",type,ifelse(subType=="","",paste("$",subType,sep="")),"<- check(cList$",type,ifelse(subType=="","",paste("$",subType,sep="")),",values)",sep="")))
  }
  else{
    ## Object initialisation: set list
    eval(parse(text=paste("cList$",type,ifelse(subType=="","",paste("$",subType,sep="")),"<- values",sep="")))
  }
  
  cList <- as.list(cList)
  return(cList)
  ##cListType <- as.list(cListType)
  ##return(cListType)
}

