## internal functions
## TODO check GenomicFeatures

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
".readGffGtf" <- function(filename=filename,annotation.type=c("exon"),ignoreWarnings=FALSE,...){
  
  ## check for the annotation type
  ## TODO we only accept one, does that make sense?
  annotation.type <- match.arg(annotation.type)

  if(!ignoreWarnings){
    if(length(list(...)) != 0){
      warning("In .readGffGtf: Ignoring extra argument(s).")
    }
  }
  
  ## sanity check
  if(!file.exists(filename)){
    stop(paste("The filename you provided:",filename,"does not exists"))
  }
  
  ## gff 3?
  if(sub("\\D+","",readLines(filename,1))!=3){
    stop(paste("Your file:",filename,"does not contain a gff header: '##gff-version 3' as first line. Is that really a gff3 file?"))
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
	all.annotation <- .readGffGtf(filename=filename,ignoreWarnings=ignoreWarnings,...)
        
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

".getGtfRange" <- function(organism=character(1),filename=filename,ignoreWarnings=FALSE,...){
	
	## read the file and do sanity checks
	all.annotation <- .readGffGtf(filename=filename,ignoreWarnings=ignoreWarnings,...)
	
	## extract the attributes
	gffAttr <- do.call(rbind,strsplit(all.annotation$gffAttributes," "))
	
	## we need 3,5, 7 and 9
	## gene 
	all.annotation$gene <- substr(gffAttr[,3],2,19)
	
	## transcript
	all.annotation$transcript <- substr(gffAttr[,5],2,19)
	
	## exon
	all.annotation$exon <- paste(all.annotation$gene,gsub("\";?","",gffAttr[,7]),sep="_")
	
	## gene name
	all.annotation$gene.name <- gsub("\";?","",gffAttr[,9])
	
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

