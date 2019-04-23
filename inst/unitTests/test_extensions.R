## unit tests for all easyRNASeq extension of other packages

### =========================
## GenomicRanges
### =========================
"test_GRangesColnames" <- function(){

  gAnnot <- RangedData(
    IRanges(
      start=c(10,30,100),
      end=c(21,53,123)),
    space=c("chr01","chr01","chr02"),
    strand=c("+","+","-"),
    transcript=c("trA1","trA2","trB"),
    gene=c("gA","gA","gB"),
    exon=c("e1","e2","e3")
  )

  exptCols <- c("transcript","gene","exon")

  ## an example of a GRangesList annotation
  grngs <- as(gAnnot,"GRanges")

  ## accessing the colnames
  checkEquals(colnames(grngs),exptCols)

  ## creating a GRangesList
  grngsList<-split(grngs,seqnames(grngs))

  ## accessing the colnames
  checkEquals(colnames(grngsList),exptCols)

}

"test_BamFileValidate" <- function(){

  tdir <- tutorialData()

  filenames <- dir(tdir,pattern="[A,C,T,G]{6}\\.bam$",full.names=TRUE)

  indexnames <- sapply(paste0(sub(".*_","",basename(filenames)),".bai"),fetchData)

  bfl <- BamFileList(filenames,index=indexnames)

  checkTrue(validate(bfl))

}