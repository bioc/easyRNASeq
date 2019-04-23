### =======================
## Internal
### =======================
## RUnit does not run within the package NAMESPACE
## as we trimline the import and export
## we need to take that into account
## hence the internal bit added to the file name

## to test the internal annotation methods
"test_internal_getAnnotation" <- function(){

  ### =======================
  ## biomaRt
  ### =======================
  if(Biobase::testBioCConnection()){
    annotParam <- AnnotParam(datasource="Dmelanogaster",type="biomaRt")
    checkTrue(easyRNASeq:::.validate(annotParam))
  } else {
    warning("No internet connection available for 'test_internal_getAnnotation' type 'BiomaRt'")
  }

  ### =======================
  ### object
  ### =======================
  grngs <- GRanges(
      ranges=IRanges(
      start=c(10,30,100),
      end=c(21,53,123)),
      seqnames=c("chr01","chr01","chr02"),
      strand=c("+","+","-"),
      transcript=c("trA1","trA2","trB"),
      gene=c("gA","gA","gB"),
      exon=c("e1","e2","e3"))

  annotParam <- AnnotParam(datasource=grngs)
  checkTrue(easyRNASeq:::.validate(annotParam))

  ### =======================
  ### rda
  ### =======================
  annotParam <- AnnotParam(fetchData("gAnnot.rda"),type="rda")
  checkTrue(easyRNASeq:::.validate(annotParam))

  ### =======================
  ### gtf
  ### =======================
  annotParam <- AnnotParam(fetchData("Drosophila_melanogaster.BDGP5.77.with-chr.gtf.gz"),type="gtf")
  checkTrue(easyRNASeq:::.validate(annotParam))

  ### =======================
  ### gff3
  ### =======================
  annotParam <- AnnotParam(datasource=fetchData("Dmel-mRNA-exon-r5.52.gff3.gz"))
  checkTrue(easyRNASeq:::.validate(annotParam))

}

### =======================
## Public
### =======================
## RUnit does not run within the package NAMESPACE
## as we trimline the import and export
## we need to take that into account
## These are public function and should not
## fail because of missing exports!

## to test the annotation methods
"test_getAnnotation_BiomaRt" <- function(){

  ## make sure we have connection
  if(Biobase::testBioCConnection()){

    ## a flybase accession
    FbAcc <- "FBgn0011656"

    ## create a biomaRt annot param
    grngs <- getAnnotation(AnnotParam(
      datasource="Dmelanogaster",
      type="biomaRt"),
                           filters="ensembl_gene_id",
                           values=FbAcc
    )

    ## check
    checkTrue(ncol(elementMetadata(grngs))==3)
    checkEquals(unique(grngs$gene),FbAcc)
  } else {
    warning("No internet connection available for 'test_getAnnotation_BiomaRt'")
  }
}

"test_getAnnotation_env" <- function(){
  expected <- GRanges(
      ranges=IRanges(
        start=c(10,30,100),
        end=c(21,53,123)),
      seqnames=c("chr01","chr01","chr02"),
      strand=c("+","+","-"),
      transcript=c("trA1","trA2","trB"),
      gene=c("gA","gA","gB"),
      exon=c("e1","e2","e3"))

  obtained <- getAnnotation(AnnotParam(
    datasource=expected))

  checkIdentical(expected,obtained)
}

"test_getAnnotation_gff3" <- function(){
  obtained <- getAnnotation(AnnotParam(
    datasource=fetchData("Dmel-mRNA-exon-r5.52.gff3.gz")))
  checkEquals(length(obtained),177816)
  checkIdentical(levels(seqnames(obtained)),
                 paste("chr",c("2L","2LHet","2R","2RHet",
                               "3L","3LHet","3R","3RHet",
                               "4","M","U","X","XHet","YHet"),sep=""))
}

"test_getAnnotation_gtf" <- function(){
  ## TODO implement me once there's a dataset

}

"test_getAnnotation_rda" <- function(){
  rngData <- getAnnotation(AnnotParam(
    datasource=fetchData("gAnnot.rda"),type="rda"))

  ## hardcoded tests
  checkTrue(ncol(elementMetadata(rngData))==3)
  checkTrue(length(rngData)==110595)
  checkIdentical(colnames(rngData),c("exon","transcript","gene"))
}
