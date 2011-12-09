### TODO for the exonCounts and other counts, we really need to pay attention to the order of the chromosomes.
### We could have that forced (i.e. having them ordered) at the rnaSeq creation time
### and have a validity action to re-order them any time the field get changed
### TODO check if the above is true or if it is already handled by IRanges and co.
### TODO think of variable width reads

### Count methods
## exons
## TODO maybe for reporting, we should simplify the list
## at least mention that some exons are duplicated
## TODO see if we can extract the aggregate to an internal f and use it for exon and island
## TODO check with alejandro and simon for exons
## TODO think if we can add a function for the user to calculate the coverage

## TODO think if we can have the genomicAnnotation re-setting somewhere else

setMethod(
          f="exonCounts",
          signature="RNAseq",
          definition=function(obj){

            genomicAnnotation(obj) <- genomicAnnotation(obj)[names(ranges(obj)) %in% names(readCoverage(obj))]
            exCounts <- .doBasicCount(obj)
            names(exCounts) <- .getName(obj,"exons")
            readCounts(obj)<-.extendCountList(readCounts(obj),exCounts,"exons",filename=fileName(obj))

            ## return
            return(obj)
          })

## features
setMethod(
          f="featureCounts",
          signature="RNAseq",
          definition=function(obj){

            genomicAnnotation(obj) <- genomicAnnotation(obj)[names(ranges(obj)) %in% names(readCoverage(obj))]
            fCounts <- .doBasicCount(obj)
            names(fCounts) <-  .getName(obj,"features")
            readCounts(obj)<-.extendCountList(readCounts(obj),
                                             fCounts,
                                             "features",
                                             filename=basename(fileName(obj)))

            ## return
            return(obj)
          })

## transcripts
setMethod(
          f="transcriptCounts",
          signature="RNAseq",
          definition=function(obj,from="exons"){

            ## check if the from is valid
            if(! from %in% c("exons","features")){
              stop("To work, the 'from' arguments should be one of 'exons', 'features'")
            }
            
            ## check if we have the exons summary already
            if(is.null(readCounts(obj,from))){
              obj <- switch(from,
                            "exons" = exonCounts(obj),
                            "features" = featureCounts(obj)
                            )
            }
            
            ## summarize these counts
            tAgg <- stats:::aggregate(readCounts(obj,from),list(transcript=.getName(obj,"transcripts")),sum)
            tCounts<-tAgg[,-1,drop=FALSE]
            rownames(tCounts)<-tAgg$transcript
            readCounts(obj)<-.extendCountList(readCounts(obj),tCounts,"transcripts",filename=fileName(obj))

            ## return
            return(obj)
          })

## genes
setMethod(
          f="geneCounts",
          signature="RNAseq",
          definition=function(obj,summarization=c("bestExons","geneModels"),...){
            
            ## check the summarization methods
            summarizations <- eval(formals("geneCounts")$summarization)
            if(!summarization %in% summarizations){
              stop(paste(
                         "The given summarization:",
                         summarization,
                         "is not part of the supported summarizations:",
                         paste(summarizations,collapse=", ")))
            }
            
            ## switch
            gCounts <- switch(EXPR=summarization,
                              "bestExons"={.bestExonSummarization(obj)},
                              "geneModels"={

                                ## check if we can at all run
                                if(length(readCoverage(obj))==0){
                                  if(length(fileName(obj))>0){
                                    stop("Cannot execute the 'geneCounts' on a multiple sample 'RNAseq' object yet. Re-run the 'easyRNASeq' function with the 'count' argument 'genes' and the 'summarization' argument 'geneModels'.")
                                  } else {
                                    stop("No coverage information available. Use either the 'easyRNASeq' or 'fetchCoverage' methods first.")
                                  }    
                                }
                                
                                ## check if we have the gene model already
                                if(nrow(geneModel(obj))==0){
                                  geneModel(obj) <- .geneModelAnnotation(genomicAnnotation(obj),...)
                                }
                                .geneModelSummarization(obj)
                              })
            readCounts(obj)<-.extendCountList(readCounts(obj),gCounts,"genes",subType=summarization,filename=fileName(obj))
            
            ## return
            return(obj)
          })

## islands
setMethod(
          f="islandCounts",
          signature="RNAseq",
          definition=function(obj,force=FALSE,...){
            
            ## check if readIsland exists, otherwise compute it
            if(is.null(readIslands(obj)$names)|force){
              obj <- findIslands(obj,...)
            }
            
            ## summarize using it
            islands=as.integer(unlist(
              aggregate(
                        readCoverage(obj)[match(names(readIslands(obj)),names(readCoverage(obj)))],
                        ranges(readIslands(obj)),
                        sum
                        )))
            
            iCounts <- ceiling(islands/readLength(obj))
            names(iCounts) <- readIslands(obj)$names
            readCounts(obj)<-.extendCountList(readCounts(obj),iCounts,"islands",filename=fileName(obj))		

            ## return
            return(obj)
          })
