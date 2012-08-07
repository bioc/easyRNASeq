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
##' @param format describes the kind of annotation provided. One of gtf or gff.
##' @param organism Organism name
##' @param gAnnot a \code{\linkS4class{RangedData}} object containing exon
##' annotations
##' @param nbCore number of CPU cores to use
##' @param filename filename that contains the annotations
##' @param \dots Additional arguments, passed to more internal functions.
##' @return A \code{\linkS4class{RangedData}} containing the loaded or
##' processed annotations.
##' @author Nicolas Delhomme
##' @keywords internal
## internal functions
## get the annot from biomaRt
".getBmRange" <- function(organism=character(1),...){
  
  ## connect
  ensembl <- useMart("ensembl")
  
  ## check the organism
  if(organism==character(1)){
    stop(paste("To use the biomaRt functionnalities, we need an organism name. Set it using the organism() function."))
  }
  
  dataset=paste(tolower(organism),"gene_ensembl",sep="_")
  
  if(! dataset %in% listDatasets(ensembl)$dataset){
    stop(paste("The organism",organism,"is not supported by the ensembl biomaRt."))
  }
  
  ## select the dataset
  ensembl=useDataset(dataset,ensembl)
  
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
  return(RangedData(
                    IRanges(
                            start=exon.annotation$exon_chrom_start,
                            end=exon.annotation$exon_chrom_end),
                    space=exon.annotation$chromosome,
                    strand=exon.annotation$strand,
                    exon=exon.annotation$ensembl_exon_id,
                    transcript=exon.annotation$ensembl_transcript_id,
                    gene=exon.annotation$ensembl_gene_id,
                    universe = organism
                    ))
  
}

## get the Annot from a gff / gtf file
".readGffGtf" <- function(filename=filename,annotation.type=c("exon"),ignoreWarnings=FALSE,format=c("gtf","gff"),...){
  
  ## TODO this will need some tidyness
  stopifnot(length(format)==1)
  
  ## check for the annotation type
  ## TODO we only accept one, does that make sense?
  annotation.type <- match.arg(annotation.type)
  
  ## sanity check
  if(!file.exists(filename)){
    stop(paste("The filename you provided:",filename,"does not exists"))
  }
  
  ## gff 3?
  if(format=="gff"){
    if(sub("\\D+","",readLines(filename,1))!=3){
      stop(paste("Your file:",filename,"does not contain a gff header: '##gff-version 3' as first line. Is that really a gff3 file?"))
    }
  }
  ## read it (genomeIntervals is the fastest)
  all.annotation <- readGff3(filename)
  
  ## keep the exon only
  ## uncomment if we have more annotation.type
  ## if(!is.null(annotation.type)){
  all.annotation<-all.annotation[all.annotation$type %in% annotation.type,]
  ##}
  if(nrow(all.annotation)==0){
    stop("The provided gff contains no annotation of type 'exon'. These are so far required.")
  }
  return(all.annotation)
}

## get the annot from a gff file
## this work for a gff from flybase, we might need another mapping
".getGffRange" <- function(organism=character(1),filename=filename,ignoreWarnings=FALSE,...){
	
	## read the file and do sanity checks
	all.annotation <- .readGffGtf(filename=filename,ignoreWarnings=ignoreWarnings,format="gff",...)
        
	## save the exon info
	all.annotation$exon <- getGffAttribute(all.annotation,"ID")
        if(all(is.na(all.annotation$exon))){
          stop("You gff file misses the ID key defining the exon ID in the gff attributes. The format should be 'gene:exon-number'.")
        }
        
	## get the gene
	all.annotation$gene <- do.call(rbind,strsplit(getGffAttribute(all.annotation,"ID"),":"))[,1]
        
        ## get the gene name
        gene.name.split <- strsplit(getGffAttribute(all.annotation,"Name"),":")
        sel<-sapply(gene.name.split,length)>2
        gene.name.split[sel] <- lapply(gene.name.split[sel],function(l){c(paste(l[-length(l)],collapse=":"),0)})
        all.annotation$gene.name <- do.call(rbind,gene.name.split)[,1]
        
	## get the transcript
	transcript <- getGffAttribute(all.annotation,"Parent")
        if(all(is.na(transcript))){
          stop("You gff file misses the Parent key defining the transcript ID in the gff attributes. For exons sharing several transcript, the transcript need to be comma separated.")
        }
	
	## duplicate the ranges
	selector<-rep(seq(along=transcript),sapply(strsplit(transcript,","),length))
	all.annotation <- all.annotation[selector,]
	
	## add the transcript
	all.annotation$transcript <- unlist(strsplit(transcript,","))
	
	## return the converted info
	exon.range <- as(all.annotation,"RangedData")
	universe(exon.range)<-organism
	
	return(exon.range)
}

".getGtfRange" <- function(organism=character(1),
                           filename=filename,
                           ignoreWarnings=FALSE,
                           fields=c("gene_id","transcript_id","exon_id","gene_name"),...){
	
	## read the file and do sanity checks
	all.annotation <- .readGffGtf(filename=filename,ignoreWarnings=ignoreWarnings,format="gtf",...)
	
	## extract the attributes
	gffAttr <- do.call(rbind,strsplit(all.annotation$gffAttributes," "))

        ## stop if the attributes we need are not present
        if(!all(fields %in% gffAttr[1,])){
          stop(paste("Your gtf file: ",filename," does not contain all the required fields: ",
                     paste(fields,collapse=", "),".",sep=""))
        }
        
        ## identify the columns we need
        sel <- match(fields,gffAttr[1,]) + 1
        
	## gene 
        ## if we have no "ENSG"
        last <- ifelse(length(grep("ENSG",gffAttr[,sel[1]]))==0,1000000L,19L)
	 
        ## remove possible annoyance
        all.annotation$gene <- gsub("\";?","",substr(gffAttr[,sel[1]],2,last))
        
	## transcript
	all.annotation$transcript <- gsub("\";?","",substr(gffAttr[,sel[2]],2,last))
	
	## exon
	all.annotation$exon <- paste(all.annotation$gene,gsub("\";?","",gffAttr[,sel[3]]),sep="_")
	
	## gene name
	all.annotation$gene.name <- gsub("\";?","",gffAttr[,sel[4]])
	
	## done
	exon.range <- as(all.annotation,"RangedData")
	universe(exon.range)<-organism
	return(exon.range)
}

### TODO we need to adapt the exon numbering depending on the strand!!!
## TODO we should take advantage of the disjoin function in the future
## Return a RangedData containing all genes model
".geneModelAnnotation" <- function(gAnnot,nbCore=1){
 
  ## can we parallelize
  if(nbCore>1){  
    ## set the number of cores
    cluster <- makePSOCKcluster(nbCore)
    
    ## get the gene ranges
    RL_list <- do.call("c",parLapply(cluster,names(gAnnot),function(chr,gAnnot){coverage(split(IRanges(start=start(gAnnot[chr]),end=end(gAnnot[chr])),gAnnot[chr]$gene))},gAnnot))
    
    ## stop the cluster
    stopCluster(cl=cluster)
  } else {
    RL_list <- do.call("c",lapply(names(gAnnot),function(chr,gAnnot){coverage(split(IRanges(start=start(gAnnot[chr]),end=end(gAnnot[chr])),gAnnot[chr]$gene))},gAnnot))
  }

  ## let's help free memory
  gc()
  
  ## TODO THIS TAKES TOO LONG
  ## get the synthetic exons
  RL<-(IRangesList(RL_list>0))
  sel<-rep(match(names(RL),gAnnot$gene),sapply(RL,length))
  return(RangedData(
                    IRanges(start=unlist(start(RL)),end=unlist(end(RL))),
                    space=gAnnot$space[sel],
                    gene=gAnnot$gene[sel],
                    strand=gAnnot$strand[sel],
                    transcript=paste(gAnnot$gene[sel],"transcript",sep="_"),
                    exon=paste(gAnnot$gene[sel],unlist(sapply(RL,function(gene){c(1:length(gene))})),sep="_")
                    ))
}

