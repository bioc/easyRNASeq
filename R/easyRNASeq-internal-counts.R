##' Internal count and summarization methods
##'
##' These are internal methods related to counting and summarizing reads
##' \itemize{
##' \item For counting reads:
##' \itemize{
##' \item .doCount A dispatcher higher level function to count and summarize reads.
##' Externalized so that it can be parallelized.
##' \item .doBasicCount A function to calculate the counts for 'exons' or 'features'
##' }
##' \item For summarizing per genes: these methods are called by the method \code{geneCounts}. Having performed the
##' \code{exonCounts} is a pre-requisite.
##' \itemize{
##' \item .bestExonSummarization Identify the exon showing the highest coverage.
##' \item .geneModelSummarization Sum the coverage values of the synthetic exons
##' constituting a gene model.
##' }
##' \item For managing the summarized read count structure:
##' \itemize{
##' \item .extendCountList extend or create the result count list of matrices
##' }
##' }
##' 
##' @aliases .doCount .doBasicCount .bestExonSummarization .geneModelSummarization
##' .extendCountList
##' @name easyRNASeq summarization internal methods
##' @rdname easyRNASeq-summarization-internal-methods
##' @param chr.map A data.frame describing the mapping of original chromosome
##' names towards wished chromosome names. See the details in
##' \code{\link[easyRNASeq:easyRNASeq]{easyRNASeq}}.
##' @param chr.sel A vector of chromosome names to subset the final results.
##' @param cList list of lists that contain count results
##' @param count The feature used to summarize the reads. One of
##' 'exons','features','genes','islands' or 'transcripts'.
##' @param filename The full path of the file to use
##' @param filter The filter to be applied when loading the data using the
##' "aln" format
##' @param format The format of the reads, one of "aln","bam". If not "bam",
##' all the types supported by the \pkg{ShortRead} package are supported too.
##' As of version 1.3.5, it defaults to bam.
##' @param gapped Is the bam file provided containing gapped alignments?
##' @param min.cov When computing read islands, the minimal coverage to take
##' into account for calling an island
##' @param min.length The minimal size an island should have to be kept
##' @param max.gap When computing read islands, the maximal gap size allowed
##' between two islands to merge them
##' @param obj An object derived from class \code{\linkS4class{RNAseq}}
##' @param plot Whether or not to plot assessment graphs.
##' @param rnaSeq An object derived from class \code{\linkS4class{RNAseq}}
##' @param summarization A character defining which method to use when
##' summarizing reads by genes. So far, only "geneModels" is available.
##' @param silent set to TRUE if you do not want messages to be printed out.
##' @param subType character string defining a sub type of counts, i.e. for the
##' gene type one of bestExon or geneModel
##' @param type \itemize{
##' \item .extendCountList: character string specifying the type of count ("exons",
##' "transcripts", "genes" or islands)
##' \item .doCount: the type of data when using the "aln" format. See the ShortRead
##' library.
##' }
##' @param validity.check Shall UCSC chromosome name convention be enforced?
##' This is only supported for a set of organisms, see
##' \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq:knownOrganisms}},
##' otherwise the argument 'chr.map' can be used to complement it.
##' @param values a named vector containing count results
##' @param \dots additional arguments. See the details in
##' \code{\link[easyRNASeq:easyRNASeq]{easyRNASeq}}.
##' @seealso
##' \itemize{
##' \item \code{\link[ShortRead:readAligned]{ShortRead:readAligned}}
##' \item \code{\linkS4class{RNAseq}}
##' \code{\link[easyRNASeq:easyRNASeq]{easyRNASeq}}.
##' }
##' 
##' @return \itemize{
##' \item .doCount: a list containing \itemize{
##' \item counts: the summarized counts as a matrix of dimension number of genes x 1
##' \item size: the library size
##' }
##' \item .doBasicCount: a vector containing read counts.
##' \item .bestExonSummarization: a vector containing summarized counts.
##' \item .geneModelSummarization: a vector containing summarized counts.
##' \item .extendCountList: a named list of matrices. The names are
##' according to the counting/summarization already performed.
##' }
##' @author Nicolas Delhomme
##' @keywords internal
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

".doCount" <- function(filename,
                       rnaSeq=rnaSeq,
                       format=c("bam","aln"),
                       filter=srFilter(),
                       count=c('exons','features','genes','islands','transcripts'),
                       type="SolexaExport",
                       chr.map=data.frame(),
                       chr.sel=c(),
                       validity.check=TRUE,
                       summarization=c("bestExons","geneModels"),
                       max.gap=integer(1),
                       min.cov=1L,
                       min.length=integer(1),
                       plot=TRUE,
                       gapped=FALSE,
                       silent=FALSE,...){

  ## load the library
  ## since we start fresh R instances
  library(easyRNASeq)

  ## some validity check
  ## useful only if I (the developer) misbehave.
  stopifnot(length(format)==1)
  .checkArguments("easyRNASeq","format",format)
  stopifnot(length(count)==1)
  .checkArguments("easyRNASeq","count",count)
  if(count=="genes"){
    stopifnot(length(summarization)==1)
  }
  if(length(summarization)==1){
    stopifnot(count=="genes")
    .checkArguments("easyRNASeq","summarization",summarization)
  }
  
  ## report
  if(!silent){
    .catn(paste("Processing",basename(filename)))
  }
  ## Fetch coverage
  rnaSeq <- fetchCoverage(rnaSeq,
                          format=format,
                          filename=filename,
                          filter=filter,type=type,
                          chr.sel=chr.sel,
                          validity.check=validity.check,
                          chr.map=chr.map,
                          gapped=gapped,...)
  
  ## emergency stop
  if(length(intersect(names(readCoverage(rnaSeq)),names(genomicAnnotation(rnaSeq))))==0
     & organismName(rnaSeq) == character(1)
     & validity.check == FALSE){
    stop(paste("Emergency stop.",
               "The chromosome names in your bam file do not match those in your annotation.",
               "You might solve that issue by providing a value to the 'organism' parameter and",
               "making sure that the 'validity.check' is set to 'TRUE'.",
               "Or you can select 'custom' as an organims and use the 'chr.map' argument to define",
               "the conversion to be applied to the chromosome names",sep="\n"))
  }
  
  ## Do count
  rnaSeq <- switch(count,
                "exons"=exonCounts(rnaSeq),
                "features"=featureCounts(rnaSeq),
                ## no need for the nbCore here, the gene model was already done
                "genes"=geneCounts(rnaSeq,summarization),
                "transcripts"=transcriptCounts(rnaSeq),
                "islands"=islandCounts(rnaSeq,max.gap=max.gap,
                  min.cov=min.cov,min.length=min.length,plot=plot)
                )
  
  return(list(counts=readCounts(rnaSeq,count,summarization),size=librarySize(rnaSeq)))
}


