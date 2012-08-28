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
            .catn(class(x), "annotation and data for", organismName(x), "containing:\n")
            ## annotations
            .catn("1) Annotations:\n")
            
            
            ## chromosomes
            if(verbose){
              .catn("\ta) Chromosome sizes:\n")
              if(length(chrSize(x)>3)){
                sToDisp <- unlist(chrSize(x))[c(1,2,length(chrSize(x)))]
                sToDisp[4]<-sToDisp[3]
                names(sToDisp)[4]<-names(sToDisp)[3]
                sToDisp[3]<-"..."
                names(sToDisp[3])<-"..."
                print(sToDisp)
              } else {
                .catn(unlist(chrSize(x)))
              }
            } else {
              .catn("\ta)",length(chrSize(x)),"chromosomes\n")
            }
            
            ## genomic Annotations (inc. geneModel and readIslands)
            if(verbose){
              .catn("\n\tb) Genomic annotations:\n")
              print(genomicAnnotation(x))
              .catn("\n\tc) Gene models\n")
              print(geneModel(x))
              .catn("\n\td) Read islands\n")
              print(readIslands(x))
            } else {
              .catn("\n\tb)",length(genomicAnnotation(x)),"genomic annotations\n")
              .catn("\n\tc)",length(geneModel(x)),"gene models\n")
              .catn("\n\td)",length(readIslands(x)),"read islands\n")
            }
            
            ## data
            .catn("\n\n2) Data:")
            if(length(readLength(x))>1){
              .catn("\tcoming from reads of length spanning the range: ",min(readLength(x)),"-",max(readLength(x)))
            } else {
              .catn("\tcoming from reads of length:",readLength(x))
            }
            if(verbose){
              .catn("\tresulting in the coverage:\n")
              .catn(show(readCoverage(x)))
            } else {
              if(length(chrSize(x))>0){
                if(length(readCoverage(x))>0){
                  sel<-match(names(readCoverage(x)),names(chrSize(x)))
                  .catn(
                      paste("\thaving an average ",
                            signif(
                                   mean(
                                          ## this throw an integer overflow
                                          ##                                          sum(readCoverage(x))/ unlist(chrSize(x)[sel])
                                          sapply(readCoverage(x),function(rC){sum(as.numeric(runLength(rC) * runValue(rC)))}
                                                 )/ unlist(chrSize(x)[sel])
                                          
                                          
                                          ),
                                   digits=2),
                            "X coverage.\n",
                            sep=""))
                }
              }
            }
            
            ## results
            if(length(readCounts(x))>0){
              .catn("\n3) Results:")
              i<-j<-1
              section <- c("a","b","c","d")
              subsection<-c(1,2)
              for(count.name in names(readCounts(x))){
                switch(EXPR=count.name,
                       "exons"={
                         .catn("\n\t",section[i],") exon summarization",sep="")
                         .catn(str(readCounts(x)$exons))
                       },
                       "features"={
                         .catn("\n\t",section[i],") feature summarization",sep="")
                         .catn(str(readCounts(x)$features))
                       },
                       "island"={
                         .catn("\n\t",section[i],") island summarization",sep="")
                         .catn(str(readCounts(x)$islands))
                       },
                       "transcripts"={
                         .catn("\n\t",section[i],") transcript summarization",sep="")
                         .catn(str(readCounts(x)$transcripts))
                       },
                       "genes"={
                         .catn("\n\t",section[i],") gene summarization",sep="")
                         for(gene.name in names(readCounts(x)$genes)){
                           switch(EXPR=gene.name,
                                  "bestExons"={
                                    .catn("\n\t\t",subsection[j],") best exon summarization",sep="")
                                    .catn(str(readCounts(x)$genes$bestExon))
                                  },
                                  "geneModels"={
                                    .catn("\n\t\t",subsection[j],") gene model summarization",sep="")
                                    .catn(str(readCounts(x)$genes$geneModel))
                                  })
                           j<-j+1
                         }
                       })
                i<-i+1
              }
            }            
          })

