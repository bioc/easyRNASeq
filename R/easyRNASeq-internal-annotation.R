##' Internal easyRNASeq annotation methods
##' 
##' These are internal methods used to retrieve annotations tabularll
##' \code{.getBmRange}Use \code{\link[biomaRt:useMart]{biomaRt}} to get exon
##' annotations.  \code{.getGffRange}Use
##' \code{\linkS4class{Genome_intervals_stranded}} to get annotation from a gff
##' file.  \code{.getGtfRange}Use
##' \code{\linkS4class{Genome_intervals_stranded}} to get annotation from a gtf
##' file.  \code{.geneModelAnnotation}Use the provided exon annotation to
##' define gene models.  \code{.readGffGtf}Use
##' \code{\linkS4class{Genome_intervals_stranded}} to get annotation from a gff
##' or gtf file. It is called from \code{getGffRange} and \code{getGtfRange}.
##' 
##' To use multicore machines more efficiently, the default parallel package
##' will be used to parallelize the processing.
##' 
##' @aliases .getBmRange .getGffRange .getGtfRange .geneModelAnnotation
##' .readGffGtf
##' @name easyRNASeq annotation internal methods
##' @rdname easyRNASeq-annotation-internal-methods
##' @param annotation.type describes the kind of annotation to keep the
##' information from in a gtf or gff file. If set to NULL all the annotations
##' are returned.
##' @param fields added a parameter that allows defining the fields parsed from
##' a gtf file. Still internal, but could easily be externalized.
##' @param filename filename that contains the annotations
##' @param format describes the kind of annotation provided. One of gtf or gff.
##' @param gAnnot a \code{\linkS4class{RangedData}} object containing exon
##' annotations
##' @param nbCore number of CPU cores to use
##' @param obj an \code{\linkS4class{AnnotParam}} object containing the 
##' necessary retrieval information (datasource and type)
##' @param ... Additional arguments, passed to more internal functions.
##' @return A \code{\linkS4class{RangedData}} containing the loaded or
##' processed annotations.
##' @author Nicolas Delhomme
##' @keywords internal

## internal functions
### =======================
## get the annot from biomaRt
### =======================
".getBmRange" <- function(obj,...){
  
  ## for the developer
  stopifnot(is(obj,"AnnotParam"))
  
  ## connect
  ensembl <- useMart(
    biomart="ensembl",
    dataset=paste(tolower(datasource(obj)),"gene_ensembl",sep="_"))
  
  ## query
  ## we extract the necessary valid arguments from ...
  exon.annotation<-eval(parse(text=paste('getBM(
                         c("ensembl_gene_id",
                           "strand",
                           "ensembl_transcript_id",
                           "chromosome_name",
                           "ensembl_exon_id",
                           "exon_chrom_start",
                           "exon_chrom_end"),
                         mart=ensembl',.getArguments("getBM",...),")",sep="")))
  
  ## convert and return
  return(GRanges(
    seqnames=exon.annotation$chromosome,
    ranges=IRanges(
      start=exon.annotation$exon_chrom_start,
      end=exon.annotation$exon_chrom_end),
    strand=exon.annotation$strand,
    exon=exon.annotation$ensembl_exon_id,
    transcript=exon.annotation$ensembl_transcript_id,
    gene=exon.annotation$ensembl_gene_id
  ))  
}

### =======================
## get the Annot from a gff / gtf file
### =======================
".readGffGtf" <- function(obj){
  
  ## read it (genomeIntervals is the fastest)
  all.annotation <- readGff3(datasource(obj))
  
  ## keep the annotation.type matching the annotation
  return(all.annotation[all.annotation$type %in% ANNOTATION.TYPE,])  
}

### =======================
## get the annot from a gff file
### =======================
".getGffRange" <- function(obj){
	
  ## read the file and do sanity checks
  all.annotation <- .readGffGtf(obj)
  
  ## save the exon info
  exon.sel <- all.annotation$type %in% ANNOTATION.TYPE["exon"]
  mRNA.sel <- all.annotation$type %in% ANNOTATION.TYPE["mRNA"]

  exons <- all.annotation[exon.sel]
  exons$exon <- getGffAttribute(all.annotation[exon.sel],"ID")
  
  ## get the gene ID
  exons$gene <- getGffAttribute(all.annotation[mRNA.sel],"Parent")[match(
    sapply(strsplit(getGffAttribute(all.annotation[exon.sel],"Parent"),","),"[",1),
    getGffAttribute(all.annotation[mRNA.sel],"ID"))]    
  
  ## get the transcript
  transcripts <- getGffAttribute(all.annotation[exon.sel],"Parent")
  
  ## duplicate the ranges
  selector<-rep(seq(along=transcripts),sapply(strsplit(transcripts,","),length))
  exons <- exons[selector,]
  
  ## add the transcript
  exons$transcript <- unlist(strsplit(transcripts,","))
  
  ## return the converted info
  return(as(exons,"GRanges"))
}

### =======================
## get the annot from a gtf file
### =======================
".getGtfRange" <- function(obj){
  
  ## read the file and do sanity checks
  all.annotation <- .readGffGtf(obj)

  ## subset for exons
  all.annotation[all.annotation$type %in% ANNOTATION.TYPE["exon"],]
  
  ## extract the attributes
  gffAttr <- do.call(rbind,strsplit(all.annotation$gffAttributes," |;"))
  
  ## identify the columns we need
  sel <- match(GTF.FIELDS,gffAttr[1,]) + 1
  
  ## gene 
  ## if we have no "ENSG"
  last <- ifelse(length(grep("ENSG",gffAttr[,sel[1]]))==0,1000000L,19L)
  
  ## remove possible annoyance
  all.annotation$gene <- gsub(" |\"|;","",substr(gffAttr[,sel[1]],1,last))
  
  ## transcript
  all.annotation$transcript <- gsub(" |\"|;","",substr(gffAttr[,sel[2]],1,last))
  
  ## exon
  ## create the exon number if they are missing
  if(is.na(sel[3])){
    exonNum <- lapply(runLength(Rle(all.annotation$gene)),":",1)
    sel <- strand(all.annotation[match(unique(all.annotation$gene),all.annotation$gene)]) == "+"
    exonNum[sel] <- lapply(exonNum[sel],rev)
    all.annotation$exon <- paste(all.annotation$gene,unlist(exonNum),sep="_")
  } else {
    all.annotation$exon <- paste(all.annotation$genes,gsub(" |\"|;","",gffAttr[,sel[3]]),sep="_")
  }
  
  ## gene name
  ## we can only have one NA: gene_name, if so, get the gene_id instead
  if(is.na(sel[4])){
    all.annotation$gene.name <- all.annotation$gene
  } else {
    all.annotation$genes.name <- gsub(" |\"|;","",gffAttr[,sel[4]])
  }
  
  ## done
  return(as(all.annotation,"GRanges"))
}

### =======================
### TODO we need to adapt the exon numbering depending on the strand!!!
## TODO we should take advantage of the disjoin function in the future
## Return a RangedData containing all genes model
".geneModelAnnotation" <- function(gAnnot,nbCore=1){
 
  ## can we parallelize
  if(nbCore>1){  
    ## set the number of cores
    cluster <- makePSOCKcluster(nbCore)
    
    ## get the gene ranges
    RL_list <- do.call("c",parLapply(cluster,names(gAnnot),function(chr,gAnnot){coverage(split(IRanges(start=start(gAnnot[chr]),end=end(gAnnot[chr])),gAnnot[chr]$genes))},gAnnot))
    
    ## stop the cluster
    stopCluster(cl=cluster)
  } else {
    RL_list <- do.call("c",lapply(names(gAnnot),function(chr,gAnnot){coverage(split(IRanges(start=start(gAnnot[chr]),end=end(gAnnot[chr])),gAnnot[chr]$genes))},gAnnot))
  }

  ## let's help free memory
  gc()
  
  ## TODO THIS TAKES TOO LONG
  ## get the synthetic exons
  RL<-(IRangesList(RL_list>0))
  sel<-rep(match(names(RL),gAnnot$genes),sapply(RL,length))
  return(RangedData(
                    IRanges(start=unlist(start(RL)),end=unlist(end(RL))),
                    space=gAnnot$space[sel],
                    gene=gAnnot$gene[sel],
                    strand=gAnnot$strand[sel],
                    transcript=paste(gAnnot$gene[sel],"transcript",sep="_"),
                    exon=paste(gAnnot$gene[sel],unlist(sapply(RL,function(gene){c(1:length(gene))})),sep="_")
                    ))
}

