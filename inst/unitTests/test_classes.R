## test the classes defaults
### =========================
## AnnotParam related
### =========================
"test_AnnotParamCharacter" <- function(){

  obj <- new("AnnotParamCharacter")

  ## the object fields
  checkTrue(is(obj@datasource,"character"))
  checkTrue(is(obj@type,"character"))

  ## the object defaults
  checkTrue(length(obj@datasource) == 0)

  ## with different parameters
  obj <- new("AnnotParamCharacter",datasource="datasource",type="type")
  checkTrue(is(obj@datasource,"character"))
  checkTrue(is(obj@type,"character"))
  checkTrue(length(obj@datasource) == 1)
  checkTrue(length(obj@type) == 1)
  checkEquals(obj@type, "type")
}

"test_AnnotParamObject" <- function(){

  obj <- new("AnnotParamObject")

  ## the object fields
  checkTrue(is(obj@datasource,"Vector"))
  checkTrue(is(obj@datasource,"GRangesList"))

  ## the object defaults
  checkTrue(length(obj@datasource) == 0)

}

### =========================
## BamParam
### =========================
"test_BamParam" <- function(){
  obj <- new("BamParam")

  ## the object fields
  checkTrue(is(obj@paired,"logical"))
  checkTrue(is(obj@stranded,"logical"))
  checkTrue(is(obj@strandProtocol,"character"))
  checkTrue(is(obj@yieldSize,"integer"))

  ## the object defaults
  checkTrue(obj@paired)
  checkTrue(!obj@stranded)
  checkTrue(obj@yieldSize==1e6L)
}

### =========================
## RnaSeqParam
### =========================

"test_RnaSeqParam" <- function(){
  obj <- new("RnaSeqParam")

  ## the object fields
  checkTrue(is(obj@countBy,"character"))
  checkTrue(is(obj@annotParam,"AnnotParam"))
  checkTrue(is(obj@bamParam,"BamParam"))
  checkTrue(is(obj@precision,"character"))

  ## the object defaults
  checkTrue(obj@countBy=="transcripts")
  checkTrue(obj@precision=="read")
}

### =========================
## TODO add the RNAseq class (or not, we may just deprecated it soon)
### =========================
