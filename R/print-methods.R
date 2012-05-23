##' Method to print a RNAseq object
##' 
##' Print information about a \code{\linkS4class{RNAseq}} object.
##' 
##' 
##' @name print methods
##' @rdname print-methods
##' @aliases print
##' @docType methods
##' @param rnaSeq An object derived from class \code{\linkS4class{RNAseq}}
##' @param verbose A logical to have a verbose or not output. Default to FALSE
##' @param \dots Additional arguments, currently unused.
##' @return Print information about a \code{\linkS4class{RNAseq}} object.
##' @author Nicolas Delhomme
##' @keywords methods
## TODO add the islands
setMethod(
          f="print",
          signature="RNAseq",
          definition=function(x,verbose=FALSE,...){
            cat(class(x), "annotation and data for", organismName(x), "containing:\n\n")
            ## annotations
            cat("1) Annotations:\n\n")
            
            
            ## chromosomes
            if(verbose){
              cat("\ta) Chromosome sizes:\n\n")
              if(length(chrSize(x)>3)){
                sToDisp <- unlist(chrSize(x))[c(1,2,length(chrSize(x)))]
                sToDisp[4]<-sToDisp[3]
                names(sToDisp)[4]<-names(sToDisp)[3]
                sToDisp[3]<-"..."
                names(sToDisp[3])<-"..."
                print(sToDisp)
              } else {
                cat(unlist(chrSize(x)),"\n")
              }
            } else {
              cat("\ta)",length(chrSize(x)),"chromosomes\n\n")
            }
            
            ## genomic Annotations (inc. geneModel and readIslands)
            if(verbose){
              cat("\n\tb) Genomic annotations:\n\n")
              print(genomicAnnotation(x))
              cat("\n\tc) Gene models\n\n")
              print(geneModel(x))
              cat("\n\td) Read islands\n\n")
              print(readIslands(x))
            } else {
              cat("\n\tb)",length(genomicAnnotation(x)),"genomic annotations\n\n")
              cat("\n\tc)",length(geneModel(x)),"gene models\n\n")
              cat("\n\td)",length(readIslands(x)),"read islands\n\n")
            }
            
            ## data
            cat("\n\n2) Data:\n")
            cat("\tcoming from reads of length:",readLength(x),"\n")
            if(verbose){
              cat("\tresulting in the coverage:\n\n")
              cat(show(readCoverage(x)),"\n")
            } else {
              if(length(chrSize(x))>0){
                if(length(readCoverage(x))>0){
                  sel<-match(names(readCoverage(x)),names(chrSize(x)))
                  cat(
                      paste("\thaving an average ",
                            signif(
                                   mean(
                                          ## this throw an integer overflow
                                          ##                                          sum(readCoverage(x))/ unlist(chrSize(x)[sel])
                                          sapply(readCoverage(x),function(rC){sum(as.numeric(runLength(rC) * runValue(rC)))}
                                                 )/ unlist(chrSize(x)[sel])
                                          
                                          
                                          ),
                                   digits=2),
                            "X coverage.\n\n",
                            sep=""))
                }
              }
            }
            
            ## results
            if(length(readCounts(x))>0){
              cat("\n3) Results:\n")
              i<-j<-1
              section <- c("a","b","c","d")
              subsection<-c(1,2)
              for(count.name in names(readCounts(x))){
                switch(EXPR=count.name,
                       "exons"={
                         cat("\n\t",section[i],") exon summarization\n",sep="")
                         cat(str(readCounts(x)$exons),"\n")
                       },
                       "features"={
                         cat("\n\t",section[i],") feature summarization\n",sep="")
                         cat(str(readCounts(x)$features),"\n")
                       },
                       "island"={
                         cat("\n\t",section[i],") island summarization\n",sep="")
                         cat(str(readCounts(x)$islands),"\n")
                       },
                       "transcripts"={
                         cat("\n\t",section[i],") transcript summarization\n",sep="")
                         cat(str(readCounts(x)$transcripts),"\n")
                       },
                       "genes"={
                         cat("\n\t",section[i],") gene summarization\n",sep="")
                         for(gene.name in names(readCounts(x)$genes)){
                           switch(EXPR=gene.name,
                                  "bestExons"={
                                    cat("\n\t\t",subsection[j],") best exon summarization\n",sep="")
                                    cat(str(readCounts(x)$genes$bestExon),"\n")
                                  },
                                  "geneModels"={
                                    cat("\n\t\t",subsection[j],") gene model summarization\n",sep="")
                                    cat(str(readCounts(x)$genes$geneModel),"\n")
                                  })
                           j<-j+1
                         }
                       })
                i<-i+1
              }
            }            
          })

