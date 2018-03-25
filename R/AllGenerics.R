###==========================
## RNAseq Class
###==========================
## accessors
## 1. getters
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
           def=function(obj,count=c("exons","features","genes","islands","transcripts"),
             summarization=c("bestExons","geneModels"),unique=FALSE){
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

###==========================
## 2. setters
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

###==========================
## pre-processing methods
###==========================
## ShortRead
##' @exportMethod demultiplex
setGeneric(name="demultiplex",
           def=function(obj,
                        barcodes=c(),
                        barcodes.qty=12,
                        barcode.length=6,
                        edition.dist=2,
                        type=c("independant","within"),
                        index.only=FALSE,
                        mc.cores=1L){
             standardGeneric("demultiplex")
           })
##' @exportMethod barcodePlot
setGeneric(name="barcodePlot",
           def=function(obj,
                        barcodes=c(),
                        type=c("independant","within"),
                        barcode.length=6,
                        show.barcode=20,
                        ...){
             standardGeneric("barcodePlot")
           })

###==========================
## Annotations

##' @exportMethod getAnnotation
setGeneric(
  name="getAnnotation",
  def=function(obj,...){
  standardGeneric("getAnnotation")
})

##' @exportMethod createSyntheticTranscripts
setGeneric(
  name="createSyntheticTranscripts",
  def=function(obj,
               features = c("mRNA", "miRNA", "tRNA", "transcript"),
               verbose = TRUE, ...){
    standardGeneric("createSyntheticTranscripts")
})

setGeneric(
           name="findIslands",
           def=function(obj,max.gap=integer(1),min.cov=1L,min.length=integer(1),plot=TRUE,...){
             standardGeneric("findIslands")
           })

###==========================
## count methods
###==========================
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

###==========================
## summary methods
###==========================
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

###==========================
## coverage methods
###==========================
##' @exportMethod fetchCoverage
setGeneric(
           name="fetchCoverage",
           def=function(obj,format=c("aln","bam"),
             filename=character(1),filter=srFilter(),
             type="SolexaExport",chr.sel=c(),
             validity.check=TRUE,chr.map=data.frame(),
             ignoreWarnings=FALSE,gapped=TRUE,
               bp.coverage=FALSE,...){
             standardGeneric("fetchCoverage")
           })

###==========================
## easy processing
###==========================
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
             outputFormat=c("matrix","SummarizedExperiment","DESeq","edgeR","RNAseq"),
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

###==========================
## edgeR extension
###==========================
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

###==========================
## edgeR & DESeq extension
###==========================
##' @exportMethod plotDispersionEstimates
setGeneric(
           name="plotDispersionEstimates",
           def=function(obj,cond=NULL,log="xy",...){
             standardGeneric("plotDispersionEstimates")
           })

##' @exportMethod multivariateConditions
setGeneric(
           name="multivariateConditions",
           def=function(obj){
             standardGeneric("multivariateConditions")
           })

##' @exportMethod plotDispLSD
setGeneric(
    name="plotDispLSD",
    def=function(obj, name = NULL, ymin,
        linecol = "#00000080", xlab = "mean of normalized counts",
        ylab = "dispersion", log = "xy", cex = 0.45,...){
        standardGeneric("plotDispLSD")
    })

###==========================
## parallel extension
###==========================
##' @exportMethod parallelize
setGeneric(name="parallelize",
           def=function(obj=list(),
             fun=NULL,
             nnodes=integer(1),...){
             standardGeneric("parallelize")})

###==========================
## BamFileList
###==========================
##' @exportMethod getBamFileList
setGeneric(name="getBamFileList",
           def=function(filenames=character(0)){
             standardGeneric("getBamFileList")
           })

##' @exportMethod validate
setGeneric(
  name="validate",
  def=function(obj,header=TRUE,cross.validation=TRUE){
    standardGeneric("validate")
  })

###==========================
## GenomicRanges extension
###==========================
##' @exportMethod colnames

##' @exportMethod unsafeAppend
setGeneric(
  name="unsafeAppend",
  def=function(obj1,obj2){
    standardGeneric("unsafeAppend")
  })

###==========================
### RnaSeqParam
###==========================
##' @exportMethod RnaSeqParam
setGeneric(name="RnaSeqParam",
           def=function(
             annotParam=AnnotParam(),
             bamParam=BamParam(),
             countBy=c("exons","features","genes","transcripts"),
             precision=c("read","bp")){
             standardGeneric("RnaSeqParam")
           })

##' @exportMethod annotParam
setGeneric(name="annotParam",
           def=function(object){
             standardGeneric("annotParam")
           })

##' @exportMethod bamParam
setGeneric(name="bamParam",
           def=function(object){
             standardGeneric("bamParam")
           })

##' @exportMethod countBy
setGeneric(name="countBy",
           def=function(object){
             standardGeneric("countBy")
           })

##' @exportMethod precision
setGeneric(name="precision",
           def=function(object){
             standardGeneric("precision")
           })

###==========================
### BamParam
###==========================
##' @exportMethod BamParam
setGeneric(name="BamParam",
           def=function(
                        paired=TRUE,
                        stranded=FALSE,
                        strandProtocol=c("reverse","forward"),
                        yieldSize=1e6L){
             standardGeneric("BamParam")
           })

##' @exportMethod paired
setGeneric(name="paired",
           def=function(object){
             standardGeneric("paired")
           })

##' @exportMethod stranded
setGeneric(name="stranded",
           def=function(object){
             standardGeneric("stranded")
           })

##' @exportMethod strandProtocol
setGeneric(name="strandProtocol",
           def=function(object){
               standardGeneric("strandProtocol")
           })

## imported but need to be exported
##' @exportMethod yieldSize

###==========================
### AnnotParam
###==========================
##' @exportMethod AnnotParam
setGeneric(name="AnnotParam",
           def=function(
             datasource=character(0),
             ...){
             standardGeneric("AnnotParam")
           })

##' @exportMethod datasource
setGeneric(name="datasource",
           def=function(object){
             standardGeneric("datasource")
           })

## imported but need to be exported
##' @exportMethod type

###==========================
## simpleRNASeq
###==========================
##' @exportMethod simpleRNASeq
setGeneric(name="simpleRNASeq",
          def=function(
          bamFiles=BamFileList(),
          param=RnaSeqParam(),
          nnodes=1,
          verbose=TRUE,
          override=FALSE){
            standardGeneric("simpleRNASeq")
           })
