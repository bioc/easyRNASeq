## RUnit does not run within the package NAMESPACE
## as we trimline the import and export
## we need to take that into account
## hence the internal bit added to the file name

## test the accessors
### =========================
## AnnotParam Accessors
### =========================
"test_internal_AnnotParam_accessor" <- function(){

  ## because it's an internal function
  ## and uses functions that are not
  ## exported by easyRNASeq
  require(GenomicRanges)

  ## create a custom object
  ## does not make sense though to have
  ## an rda file with the default (gff3) type
  ## but that's just a test
  obj <- new("AnnotParamCharacter",datasource=fetchData("gAnnot.rda"))

  ## test the type default
  checkTrue(type(obj)=="gff3")

  ## check that the datasource exist
  checkTrue(file.exists(datasource(obj)))

  ## and now an AnnotParamObject
  obj <- new("AnnotParamObject",
             datasource=GRangesList())
  checkTrue(is.null(type(obj)))
  checkEquals(datasource(obj),GRangesList())
}

### =========================
## BamParam Accessors
### =========================
"test_BamParam_accessor" <- function(){

  ## create a custom object
  obj <- new("BamParam",yieldSize=1L,
             paired=FALSE,stranded=TRUE,
             strandProtocol="forward")

  checkTrue(yieldSize(obj)==1L)
  checkTrue(!paired(obj))
  checkTrue(stranded(obj))
  checkTrue(strandProtocol(obj)=="forward")
}

### =========================
## RnaSeqParam Accessors
### =========================
"test_RnaSeqParam_accessor" <- function(){

  ## create a custom object
  obj <- new("RnaSeqParam",bamParam=BamParam(yieldSize=1L))

  ## test the yieldSize inheritance
  checkTrue(yieldSize(obj)==1L)
  checkTrue(yieldSize(obj)!=1.2)

  ## test the paired and stranded inheritance
  ## of the default values
  checkTrue(paired(obj))
  checkTrue(!stranded(obj))

  ## check the defaults
  checkTrue(length(datasource(obj))==0)
  checkTrue(obj@precision=="read")
  checkTrue(obj@countBy=="transcripts")

}
