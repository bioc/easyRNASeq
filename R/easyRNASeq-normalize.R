##' easyRNASeq count table correction to RPKM
##' 
##' Convert a count table obtained from the easyRNASeq function into an RPKM
##' corrected count table.
##' 
##' RPKM accepts two sets of arguments:
##' \itemize{
##' \item{RNAseq,character}{ the
##' \dots{} are additional arguments to be passed to the
##' \code{\link[easyRNASeq:easyRNASeq-accessors]{readCounts}} method.}
##' \item{matrix,named vector}{normalize a count matrix by providing the feature
##' sizes (e.g. gene sizes) as a named vector where the names match the row
##' names of the count matrix and the lib sizes as a named vector where the
##' names match the column names of the count matrix.}}
##' 
##' @aliases RPKM RPKM,RNAseq,ANY,ANY,ANY-method
##' RPKM,matrix,ANY,vector,vector-method
##' @name easyRNASeq correction methods
##' @rdname easyRNASeq-correction-methods
##' @param feature.size Precise the feature (e.g. exons, genes) sizes. It
##' should be a named numeric list, named after the feature names.
##' @param from Determine the kind of coverage to use, choice limited to:
##' exons, features, transcripts, bestExons, geneModels or islands.
##' @param lib.size Precise the library size. It should be a named numeric
##' list, i.e. named after the sample names.
##' @param obj An object of class \code{\linkS4class{RNAseq}} or a
##' \code{matrix}, see details
##' @param simplify If set to TRUE, whenever a feature (exon, feature, ...) is
##' duplicated in the count table, it is only returned once.
##' @param ... additional arguments. See details
##' @return A \code{matrix} containing RPKM corrected read counts.
##' @author Nicolas Delhomme
##' @seealso \code{\link[easyRNASeq:easyRNASeq-accessors]{readCounts}}
##' @keywords methods
##' @examples
##' 
##' 
##' 	\dontrun{
##' 	## get an RNAseq object
##' 	rnaSeq <- easyRNASeq(filesDirectory=
##' 		    			system.file(
##' 					"extdata",
##' 					package="RnaSeqTutorial"),
##' 					pattern="[A,C,T,G]{6}\.bam$",
##' 				format="bam",
##' 				readLength=36L,
##' 				organism="Dmelanogaster",
##' 				chr.sizes=as.list(seqlengths(Dmelanogaster)),
##' 				annotationMethod="rda",
##' 				annotationFile=system.file(
##' 				                            "data",
##' 							    "gAnnot.rda",
##' 							    package="RnaSeqTutorial"),
##' 				count="exons",
##' 				outputFormat="RNAseq")
##' 
##' 	## get the RPKM
##' 	rpkm <- RPKM(rnaSeq,from="exons")
##' 	
##' 	## the same from a count table
##' 	count.table <- readCounts(rnaSeq,count="exons")
##' 
##' 	## get the RPKM
##' 	## verify that the feature are sorted as the count.table
##' 	all(.getName(rnaSeq,"exon") == rownames(count.table))
##' 	feature.size <- unlist(width(ranges(rnaSeq)))
##' 
##' 	## verify that the samples are ordered in the same way
##' 	all(names(librarySize(rnaSeq)) == colnames(count.table))
##' 
##' 	## get the RPKM
##' 	rpkm <- RPKM(count.table,
##' 			feature.size=feature.size,
##' 			lib.size=librarySize(rnaSeq))
##' }
##' 
setMethod(
          f="RPKM",
          signature=c("matrix","ANY","vector","vector"),
          definition=function(obj,from,lib.size=numeric(1),feature.size=integer(1),simplify=TRUE,...){
            
            ## SANITY check here
            if(length(lib.size) != ncol(obj)){
              stop("You need to provide a lib.size named vector. The names should be the colnames of your matrix. The vector size should be equal to the number of matrix columns.")
            }

            ## TODO sanity check for the rownames
            
            ## nothing to do with the "..." so far 
            ## maybe subselect?

            ## subselect
            if(length(feature.size) != nrow(obj)){
              stop("You need to provide a feature.size named vector. The names should be the rownames of your matrix. The vector size should be equal to the number of matrix rows.")
            }

            ## TODO sanity check for the colnames

            ## calc
            res <- do.call("cbind",lapply(
                                          c(1:ncol(obj)),function(i,obj,sizes,lib.size){
                                            (obj[,i] / (lib.size[i]/10^6)) / (sizes/10^3)
                                          },obj,feature.size[match(rownames(obj),names(feature.size))],
                                          lib.size[match(colnames(obj),names(lib.size))]))
            
            colnames(res) <- colnames(obj)
            rownames(res) <- rownames(obj)
            if(simplify){
              res<-res[!duplicated(rownames(res)),,drop=FALSE]
            }
            return(res)
          })

setMethod(
          f="RPKM",
          signature="RNAseq",
          definition=function(obj,
            from=c("exons","features","transcripts","bestExons","geneModels","islands"),
            lib.size=numeric(1),feature.size=integer(1),simplify=TRUE,...){
            
            ## get the possible values
            .checkArguments("RPKM","from",from)
            
            ## switch to get the values
            mCounts <- switch(from,
                              "exons"=readCounts(obj,'exons',...),
                              "features"=readCounts(obj,'features',...),
                              "transcripts"=readCounts(obj,'transcripts',...),
                              "bestExons"=readCounts(obj,'genes','bestExons',...),
                              "geneModels"=readCounts(obj,'genes','geneModels',...),
                              "islands"=readCounts(obj,'islands',...)
                              )
            
            ## get the lib sizes
            libSizes <- librarySize(obj)
            libSizes <- libSizes[match(colnames(mCounts),names(libSizes))]
            
            ## valid?
            if(is.null(mCounts) | length(mCounts)==0){
              stop(paste(
                         "The summarization by",
                         from,
                         "was not performed yet!",
                         "No counts can therefore be reported."
                         ))
            }
            
            ## if mCounts is a vector, matrix it
            if(is.vector(mCounts)){
              mCounts <- as.matrix(mCounts)
            }
            
            ## get the feature sizes
            ## as the feature can be filtered differently than the annotation, we need to pay attention to it too
            mSize <- switch(from,
                            "exons"=width(genomicAnnotation(obj))[match(rownames(mCounts),.getName(obj,"exons"))],
                            "features"=width(genomicAnnotation(obj))[match(rownames(mCounts),.getName(obj,"features"))],
                            "transcripts"= {
                              ## aggregate first
                              agg <- stats::aggregate(width(genomicAnnotation(obj)),list(transcript=.getName(obj,"transcripts")),sum)
                              ## then sort
                              agg[match(rownames(mCounts),agg[,1]),2]},
                            "bestExons"= width(genomicAnnotation(obj))[match(rownames(mCounts),.getName(obj,"exons"))],
                            "geneModels"= {
                              ## aggregate
                              agg<-stats::aggregate(width(geneModel(obj)),list(gene=geneModel(obj)$gene),sum)
                              ## sort
                              agg[match(rownames(mCounts),agg[,1]),2]},
                            ## same as that: as(by(width(geneModel(obj)),geneModel(obj)$gene,sum),"vector") but faster
                            "islands"= width(readIslands(obj))[match(rownames(mCounts),readIslands(obj))$names]
                            )
            
            ## we got a matrix, so work per column
            res <- do.call(cbind,lapply(c(1:ncol(mCounts)),function(i,mCounts,sizes,libSizes){
              (mCounts[,i] / (libSizes[i]/10^6)) / (sizes/10^3)
            },mCounts,mSize,libSizes))
            colnames(res) <- colnames(mCounts)
            ## filter for unique row names
            if(simplify){
              res<-res[!duplicated(rownames(res)),,drop=FALSE]
            }
            return(res)
          })

