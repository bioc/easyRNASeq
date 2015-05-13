## to test the method result consistency over time
"test_simpleRNASeq" <- function(){
  
  ## expected results
  columnSums <- c("ACACTG.bam" = 38381,
                  "ACTAGC.bam" = 28895,
                  "ATGGCT.bam" = 37569,
                  "TTGCGA.bam" = 41023)

  selectExons <- c("FBgn0003884:2"=1144,
                   "FBgn0003887:1"=1117,
                   "FBgn0038476:4"=850,
                   "FBgn0004867:1"=806,
                   "FBgn0037874:1"=709,
                   "FBgn0259745:6"=592)
  
  ## the data
  library("RnaSeqTutorial")
  
  ## get the BamFileList
  bamFiles <- getBamFileList(
    dir(path=system.file("extdata",
                         package="RnaSeqTutorial"),
        pattern="^[A,T].*\\.bam$",
        full.names=TRUE))
  
  ## create the AnnotParam
  annotParam <- AnnotParam(system.file(
    "extdata",
    "Dmel-mRNA-exon-r5.52.gff3",
    package="RnaSeqTutorial"))
  
  ## create the RnaSeqParam
  param <- RnaSeqParam(annotParam=annotParam)
  
  ## get a RangedSummarizedExperiment containing the counts table
  sexp <- simpleRNASeq(
    bamFiles=bamFiles,
    param=param,
    verbose=FALSE
  )
  
  ## check the overall counts
  checkIdentical(colSums(assay(sexp)),columnSums)
  
  ## check the selected exon counts
  checkIdentical(head(sort(rowSums(assay(sexp))[rowSums(assay(sexp))>0],
                           decreasing=TRUE)),
                 selectExons)
  
}
