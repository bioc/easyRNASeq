# get the example data and annotation files from GitHub
exFiles <- c("gAnnot.rda",
              "Drosophila_melanogaster.BDGP5.77.with-chr.gtf.gz",
              "Dmel-mRNA-exon-r5.52.gff3.gz",
              "ACACTG.bam","ACTAGC.bam","ATGGCT.bam","TTGCGA.bam",
              "ACACTG.bam.bai","ACTAGC.bam.bai","ATGGCT.bam.bai","TTGCGA.bam.bai")

invisible(sapply(exFiles,function(f){
     if(!file.exists(f)){
         invisible(download.file(paste0("https://github.com/UPSCb/UPSCb/raw/",
                                    "master/tutorial/easyRNASeq/",f),f))
     }
 }))

# run the tests
BiocGenerics:::testPackage("easyRNASeq")

# cleanup
file.remove(exFiles)
