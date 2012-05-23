##' Identify expressed regions de-novo
##' 
##' Process the coverage to locate regions with a minimum coverage (min.cov).
##' If regions are separated by a gap shorter than a maximum length (max.gap),
##' they are unified. Only islands longer than min.length are returned.  These
##' functions are now outdated and would need to be actualized.
##' 
##' \dots{} are for providing additional options to the
##' \code{\link[graphics:hist]{hist}} plot function.
##' 
##' @aliases findIslands
##' @name easyRNASeq island methods
##' @rdname easyRNASeq-island-methods
##' @param obj An object of class \code{RNAseq}
##' @param max.gap Maximum gap between two peaks to build an island
##' @param min.cov Minimum coverage for an island to be returned
##' @param min.length Minimum size of an island to be returned
##' @param plot If TRUE, draw plots of coverage distribution. Help the user to
##' select an appropriate value for the minimum coverage.
##' @param \dots See details
##' @return An \code{RNAseq} object with the readIsland slot set with a
##' RangedData containing the selected islands and the readCount slot
##' actualized with a list containing the count table per island.
##' @author Nicolas Delhomme
##' @keywords connection data methods
##' @examples
##' 
##' 	\dontrun{
##' 	## NOTE that this function might need to be actualized
##' 	obj <- new('RNAseq',
##' 		organismName="Dmelanogaster",
##' 		readLength=36L,
##' 		chrSize=as.list(seqlengths(Dmelanogaster))
##' 		)
##' 
##' 	obj <- fetchCoverage(
##' 			obj,
##' 			format="bam",
##'                         filename=system.file(
##' 				"extdata",
##' 				"ACACTG.bam",
##'                             	package="RnaSeqTutorial")
##' 			)
##' 
##' 	obj <- findIslands(
##' 			obj,
##' 			max.gap=10L,
##' 			min.cov=10L,
##' 			min.length=200L)
##' 	}
##' 
## find the islands
setMethod(
          f="findIslands",
          signature="RNAseq",
          definition=function(obj,max.gap=integer(1),min.cov=1L,min.length=integer(1),plot=TRUE,...){			
            
            ## check the gap
            if(max.gap<0){
              stop("There are no such thing as a negative gap... Change your max.gap parameter.")
            }
            
            ## check the read coverage
            if(length(readCoverage(obj))==0){
              stop("No coverage information was provided. Cannot compute the read islands")
            }
            
            ## check the min cov
            if(min.cov != 1L){
              if(is.na(as.numeric(min.cov))){
                stop("The min.cov argument accepts only an integer value")
              }
              if(!is.integer(min.cov)){
                min.cov<-as.integer(min.cov)
                warning(paste("The min.cov should be an integer, not a numeric value. Changing your value to:",min.cov))
              }
            } 
            
            ## check min.length
            if(min.length<1){
              stop("There are no such thing as reads length inferior as 1 base. Change your min.length value")
            }
            
            ## the endrule is necessary to keep the length of the coverage consistant with the chr length
            island <- IRangesList(sapply(runsum(readCoverage(obj),max.gap,endrule="constant"),function(x){IRanges(as.logical(x>0))}))
            
            ## island.length <-  IRangesList(sapply(island.gap,function(x){x[width(x)>=min.length]}))
            island <- island[width(island)>=min.length,]
            
            ## in case we have no island on a chromosome
            island <- island[sapply(island,length) != 0]
            
            ##islandsCount
            islandC <- as.integer(aggregate(
                                            readCoverage(obj)[match(names(island),names(readCoverage(obj)))],
                                            island,
                                            sum
                                            )
                                  /readLength(obj))
            
            ## plot the cov using a boxplot
            ## and an hist; hist(sort(islandC),breaks=seq(0,max(islandC),1),...)
            if(plot==TRUE){	
              par(mfrow=c(2,2))
              hist(islandC,breaks=seq(0,max(islandC),1),main='Frequency of coverage',xlab='Coverage',...)
              boxplot(islandC,log='y',main='Distribution of coverage',xlab=paste("n =",sum(sapply(island,length))))
              abline(h=median(islandC),lty=2,col='orange')
              abline(h=min.cov,lty=3,col='red')
              mtext(median(islandC),col = 'orange', side = 2,at=median(islandC), line=2)
            }
            
            ## create a selector for the islandC having a coverage >= min.cov                        
            selector <- islandC>min.cov
            
            ## subset island and islandC
            islandC <- islandC[selector]
            island <- island[LogicalList(split(selector,rep(1:length(island),sapply(island,length))))]
            
            ## plot the cov using a boxplot
            ## and an hist
            if(plot==TRUE){	
              hist(islandC,breaks=seq(0,max(islandC),1),main='Frequency of coverage after filtering',xlab='Coverage',...)
              boxplot(islandC,log='y',main='Distribution of coverage after filtering',xlab=paste("n =",sum(sapply(island,length))),sub=paste("coverage cutoff =",min.cov))
              abline(h=median(islandC),lty=2,col='orange')
              abline(h=min.cov,lty=3,col='red')
              mtext(median(islandC),col = 'orange', side = 2, at = median(islandC), line = 2)
            }
            
            ## name them
            names(islandC) <- sprintf('%s%s%07d',names(unlist(island)),'_is',c(1:(length(unlist(island)))))
            ## set the readIslands
            readIslands(obj) <- RangedData(ranges=island,names=names(islandC))
            
            ## set the counts
            readCounts(obj)<-.extendCountList(readCounts(obj),islandC,"islands",filename=fileName(obj))
            
            ## update the obj
            return(obj)
          })
