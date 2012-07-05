## accessors
## getters
##' @exportMethod genomicAnnotation
setGeneric(
           name="genomicAnnotation",
           def=function(obj){
             standardGeneric("genomicAnnotation")
           })
##' @exportMethod readCoverage
setGeneric(
           name="readCoverage",
           def=function(obj){
             standardGeneric("readCoverage")
           })
##' @exportMethod readLength
setGeneric(
           name="readLength",
           def=function(obj){
             standardGeneric("readLength")
           })
##' @exportMethod readCounts
setGeneric(
           name="readCounts",
           def=function(obj,count=c("exons","features","genes","islands","transcripts"),summarization=c("bestExons","geneModels"),unique=FALSE){
             standardGeneric("readCounts")
           }
           )
##' @exportMethod librarySize
setGeneric(
           name="librarySize",
           def=function(obj){
             standardGeneric("librarySize")
           })
##' @exportMethod organismName
setGeneric(
           name="organismName",
           def=function(obj){
             standardGeneric("organismName")
           })
##' @exportMethod geneModel
setGeneric(
           name="geneModel",
           def=function(obj){
             standardGeneric("geneModel")
           })

setGeneric(
           name="readIslands",
           def=function(obj){
             standardGeneric("readIslands")
           })
##' @exportMethod chrSize
setGeneric(
           name="chrSize",
           def=function(obj){
             standardGeneric("chrSize")
           })
##' @exportMethod fileName
setGeneric(
           name="fileName",
           def=function(obj){
             standardGeneric("fileName")
           })

## setters
##' @exportMethod genomicAnnotation<-
setGeneric(
           name="genomicAnnotation<-",
           def=function(obj,value){
             standardGeneric("genomicAnnotation<-")
           })
##' @exportMethod readLength<-
setGeneric(
           name="readLength<-",
           def=function(obj,value){
             standardGeneric("readLength<-")
           })
##' @exportMethod readCoverage<-
setGeneric(
           name="readCoverage<-",
           def=function(obj,value){
             standardGeneric("readCoverage<-")
           })
##' @exportMethod chrSize<-
setGeneric(
           name="chrSize<-",
           def=function(obj,value){
             standardGeneric("chrSize<-")
           })
##' @exportMethod readCounts<-
setGeneric(
           name="readCounts<-",
           def=function(obj,value){
             standardGeneric("readCounts<-")
           })
##' @exportMethod librarySize<-
setGeneric(
           name="librarySize<-",
           def=function(obj,value){
             standardGeneric("librarySize<-")
           })
##' @exportMethod organismName<-
setGeneric(
           name="organismName<-",
           def=function(obj,value){
             standardGeneric("organismName<-")
           })
##' @exportMethod geneModel<-
setGeneric(
           name="geneModel<-",
           def=function(obj,value){
             standardGeneric("geneModel<-")
           })

setGeneric(
           name="readIslands<-",
           def=function(obj,value){
             standardGeneric("readIslands<-")
           })
##' @exportMethod fileName<-
setGeneric(
           name="fileName<-",
           def=function(obj,value){
             standardGeneric("fileName<-")
           })

## pre-processing methods
## ShortRead
##' @exportMethod demultiplex
setGeneric(name="demultiplex",
           def=function(obj,barcodes=c(),barcodes.qty=12,barcode.length=6, edition.dist=2, type=c("independant","within"),index.only=FALSE){
             standardGeneric("demultiplex")
           })
##' @exportMethod barcodePlot
setGeneric(name="barcodePlot",
           def=function(obj,barcodes=c(),type=c("independant","within"),barcode.length=6,show.barcode=20,...){
             standardGeneric("barcodePlot")
           })

## Annotations
##' @exportMethod fetchAnnotation
setGeneric(
           name="fetchAnnotation",
           def=function(obj,
             method=c("biomaRt","gff","gtf"),
             filename=character(1),
             ignoreWarnings=FALSE,
             ...){
             standardGeneric("fetchAnnotation")
           })

## setGeneric(
##            name="getCoverage",
##            def=function(obj,aln){
##              standardGeneric("getCoverage")
##            })
## TODO think of more proper parameters here
setGeneric(
           name="findIslands",
           def=function(obj,max.gap=integer(1),min.cov=1L,min.length=integer(1),plot=TRUE,...){
             standardGeneric("findIslands")
           })

## count methods
##' @exportMethod exonCounts
setGeneric(
           name="exonCounts",
           def=function(obj){
             standardGeneric("exonCounts")
           })
##' @exportMethod featureCounts
setGeneric(
           name="featureCounts",
           def=function(obj){
             standardGeneric("featureCounts")
           })
##' @exportMethod transcriptCounts
setGeneric(
           name="transcriptCounts",
           def=function(obj,from="exons"){
             standardGeneric("transcriptCounts")
           })
##' @exportMethod geneCounts
setGeneric(
           name="geneCounts",
           def=function(obj,summarization=c("bestExons","geneModels"),...){
             standardGeneric("geneCounts")
           })

setGeneric(
           name="islandCounts",
           def=function(obj,force=FALSE,...){
             standardGeneric("islandCounts")
           })

## summary methods
##' @exportMethod RPKM
setGeneric(
           name="RPKM",
           def=function(obj,
             from=c("exons",
               "features",
               "transcripts",
               "bestExons",
               "geneModels",
               "islands"),
             lib.size=numeric(1),
             feature.size=numeric(1),
             simplify=TRUE,...){
             standardGeneric("RPKM")
           })
##' @exportMethod fetchCoverage
setGeneric(
           name="fetchCoverage",
           def=function(obj,format=c("aln","bam"),
             filename=character(1),filter=srFilter(),
             type="SolexaExport",chr.sel=c(),
             isUnmappedQuery=FALSE,what=c("rname","pos","qwidth"),
             validity.check=TRUE,chr.map=data.frame(),
             ignoreWarnings=FALSE,gapped=TRUE,bp.coverage=FALSE,...){
             standardGeneric("fetchCoverage")
           })

## easy processing
##' @exportMethod easyRNASeq
setGeneric(
           name="easyRNASeq",
           def=function(filesDirectory=character(1),
             organism=character(1),
             chr.sizes=c("auto"),
             readLength=integer(1),
             annotationMethod=c("biomaRt","env","gff","gtf","rda"),
             annotationFile=character(1),
             annotationObject=RangedData(),
             format=c("bam","aln"),gapped=FALSE,
             count=c('exons','features','genes','islands','transcripts'),
             outputFormat=c("DESeq","edgeR","matrix","RNAseq"),
             pattern=character(1),filenames=character(0),nbCore=1,
             filter=srFilter(),type="SolexaExport",
             chr.sel=c(),summarization=c("bestExons","geneModels"),
             normalize=FALSE,
             max.gap=integer(1),min.cov=1L,
             min.length=integer(1),plot=TRUE,
             conditions=c(),
             validity.check=TRUE,
             chr.map=data.frame(),
             ignoreWarnings=FALSE,
             silent=FALSE,...){
             standardGeneric("easyRNASeq")
           })

## known organisms
##' @exportMethod knownOrganisms
setGeneric(
           name="knownOrganisms",
           def= function(missing){
             standardGeneric("knownOrganisms")
           })

## edgeR extension
##' @exportMethod plotNormalizationFactors
setGeneric(
           name="plotNormalizationFactors",
           def=function(
             obj=DGEList(),
             cond1=character(1),
             cond2=character(1)
             ){
             standardGeneric("plotNormalizationFactors")
           })

## edgeR & DESeq extension
##' @exportMethod plotDispersionEstimates
setGeneric(
           name="plotDispersionEstimates",
           def=function(obj,...){
             standardGeneric("plotDispersionEstimates")
           })

## DESeq extension
## accessor
## setGeneric(
##            name="fitInfo",
##            def=function(obj){
##              standardGeneric("fitInfo")
##            })

##' @exportMethod multivariateConditions
setGeneric(
           name="multivariateConditions",
           def=function(obj){
             standardGeneric("multivariateConditions")
           })

## parallel extension
##' @exportMethod parallelize
setGeneric(name="parallelize",
           def=function(obj=list(),
             fun=NULL,
             nnodes=integer(1),...){
             standardGeneric("parallelize")})
