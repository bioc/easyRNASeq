## to test the IRanges extension
"test_ranges_accessor" <- function(){
  ## create a grange
  grng <- GRanges(seqnames=Rle(factor(c("chr1","chr3","chr11"))),
                  ranges=IRanges(1:3,4:6),
                  strand=Rle(factor(c("+","-","+"))),
                  exon=c("ex1","ex2","ex3"))
    
  ## create an RNAseq object
  rnaSeq <- new('RNAseq',genomicAnnotation=grng)
  
  ## make sure the range accessor works
  checkTrue(
    identical(unlist(ranges(rnaSeq)),
              unlist(split(ranges(grng),seqnames(grng)))))
}

"test_internal_names_accessor" <- function(){  
  ## create a grange
  grng <- GRanges(seqnames=Rle(factor(c("chr1","chr3","chr11"))),
                  ranges=IRanges(1:3,4:6),
                  strand=Rle(factor(c("+","-","+"))),
                  exon=c("ex1","ex2","ex3"))

  ## create an RNAseq object
  rnaSeq <- new('RNAseq',genomicAnnotation=grng)
  
  ## make sure the name accessor works
  checkTrue(all(easyRNASeq:::.getName(rnaSeq,"exons") ==  unlist(split(grng$exon,seqnames(grng)))))
}
