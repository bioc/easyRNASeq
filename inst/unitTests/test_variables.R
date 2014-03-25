"test_global_variables" <- function(){
  checkIdentical(GTF.FIELDS,
                 c("gene_id","transcript_id","exon_number","gene_name"))
  checkIdentical(ANNOTATION.TYPE,c(mRNA="mRNA",exon="exon"))
}