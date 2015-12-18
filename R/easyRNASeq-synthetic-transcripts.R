## TODO test me
## TODO vignette me

##' Methods to create synthetic transcripts
##'
##' This function create a set of synthetic transcripts from a provided
##' annotation file in "gff3" or "gtf" format. As detailed in
##' \url{http://www.epigenesys.eu/en/protocols/bio-informatics/1283-guidelines-for-rna-seq-data-analysis},
##' one major caveat of estimating gene expression using aligned RNA-Seq reads
##' is that a single read, which originated from a single mRNA molecule, might
##' sometimes align to several features (e.g. transcripts or genes) with
##' alignments of equivalent quality. This, for example, might happen as a result
##' of gene duplication and the presence of repetitive or common domains.
##' To avoid counting unique mRNA fragments multiple times, the
##' stringent approach is to keep only uniquely mapping reads - being aware of
##' potential consequences. Not only can "multiple counting" arise from a
##' biological reason, but also from technical artifacts, introduced mostly
##' by poorly formatted gff3/gtf annotation files. To avoid this, it is best
##' practice to adopt a conservative approach by collapsing all existing
##' transcripts of a single gene locus into a "synthetic" transcript containing
##' every exon of that gene. In the case of overlapping exons, the longest
##' genomic interval is kept, i.e. an artificial exon is created. This process
##' results in a flattened transcript - a gene structure with a one (gene) to
##' one (transcript) relationship.
##'
##' The \code{createSyntheticTranscripts} function implements this, taking
##' advantage of the hierarchical structure of the gff3/gtf file. Exon
##' features are related to their transcript (parent), which themselves derives
##' from their gene parents. Using this relationship, exons are combined per gene
##' into a flattened transcript structure. Note that this might not avoid multiple
##' counting if genes overlap on opposing strands. There, only strand specific
##' sequencing data has the power to disentangle these situations.
##'
##' As gff3/gtf file can contain a large number of feature types, the
##' \code{createSyntheticTranscripts} currently only supports: \emph{mRNA},
##' \emph{miRNA}, \emph{tRNA} and \emph{transcript}. Please contact me if you
##' need additional features to be considered. Note however, that I will only
##' add features that are part of the \url{sequenceontology.org} SOFA
##' (SO_Feature_Annotation) ontology.
##'
##' @aliases createSyntheticTranscripts
##' createSyntheticTranscripts,AnnotParamCharacter-method
##' createSyntheticTranscripts,character-method
##' @rdname easyRNASeq-synthetic-transcripts
##' @param obj a \code{\linkS4class{AnnotParamCharacter}} object or the
##' annotation filename as a \code{character} string
##' @param features one or more of 'mRNA', 'miRNA', 'tRNA', 'transcript'
##' @param ... If \code{obj} is a character string, \code{input} and
##' \code{output} - see below
##' @param input the type of input, one of 'gff3' or 'gtf'
##' @param output the output type, one of 'Genome_intervals' or 'GRanges'
##' @param verbose increase the verbosity (default TRUE)
##' @return
##' Depending on the \code{obj} class.
##' \itemize{
##'   \item \code{AnnotParamCharacter}: a \code{AnnotParamObject} object
##'   \item a \code{character} filename: depending on the selected \code{output}
##'   value, a \code{\link[genomeIntervals:Genome_intervals-class]{Genome_intervals}}
##'   or a \code{\linkS4class{GRanges}} object.
##' }
##' @author Nicolas Delhomme
##' @seealso
##' \itemize{
##' \item{For the input:
##' \itemize{
##' \item \code{\linkS4class{AnnotParam}}
##' }}
##' \item{For the output:
##' \itemize{
##' \item \code{\linkS4class{AnnotParam}}
##' \item \code{\link[genomeIntervals:Genome_intervals-class]{Genome_intervals}}
##' \item \code{\linkS4class{GRanges}}
##' }}}
##' @keywords methods
##' @examples
##'
##'   \dontrun{
##'   ## the data
##'   library("RnaSeqTutorial")
##'
##'   ## get the example file
##'   library(curl)
##'   curl_download(paste0("https://microasp.upsc.se/root/upscb-public/raw/",
##'   "master/tutorial/easyRNASeq/Drosophila_melanogaster.BDGP5.77.with-chr.gtf.gz"),
##'              "Drosophila_melanogaster.BDGP5.77.with-chr.gtf.gz")
##'
##'   ## create the AnnotParam
##'   annotParam <- AnnotParam(
##'     datasource="Drosophila_melanogaster.BDGP5.77.with-chr.gtf.gz",
##'     type="gtf")
##'
##'   ## create the synthetic transcripts
##'   annotParam <- createSyntheticTranscripts(annotParam,verbose=FALSE)
##'
##'  }
##'
## TODO change the example!
##   annotParam <- AnnotParam(system.file(
##                    "extdata",
##                    "Dmel-mRNA-exon-r5.52.gff3",
##                    package="RnaSeqTutorial"))
setMethod(f = "createSyntheticTranscripts",
          signature = "AnnotParamCharacter",
          definition = function(obj,
                                features = c("mRNA", "miRNA", "tRNA", "transcript"),
                                verbose = TRUE) {

            # Check for the input type
            supported.types <- c("gff3","gtf")
            if(!type(obj) %in% supported.types){
              stop(sprintf("Only the %s types are supported by this function",
                   paste(supported.types,collapse=", ")))
            }

            # Validate the object
            tryCatch(.validate(obj,verbose=verbose),
                     error=function(e){
                        warning(paste("Your gtf file is not comprehensive;",
                                      "possibly containing only 'exon' features"))
                        message("Validating under lenient criteria")
                        .validate(obj,verbose=verbose,lenient=TRUE)
                     })

            # Create the synth. trx.
            return(
              AnnotParam(datasource =
                           createSyntheticTranscripts(datasource(obj),
                                                      features=features,
                                                      output="GRanges",
                                                      input=type(obj),
                                                      verbose=verbose)))
})

##' @rdname easyRNASeq-synthetic-transcripts
setMethod(f = "createSyntheticTranscripts",
          signature = "character",
            definition = function(obj,
                                  features = c("mRNA", "miRNA", "tRNA", "transcript"),
                                  verbose = TRUE,
                                  output = c("Genome_intervals","GRanges"),
                                  input = c("gff3","gtf")
                                  ) {

  # first check
  stopifnot(file.exists(obj))

  # get the values
  input <- match.arg(input)
  features <- match.arg(features, several.ok = TRUE)
  output <- match.arg(output)

  # define some global variables
  relation <- switch(input,
                     "gff3"=list(ID="ID",Parent="Parent"),
                     "gtf"=list(ID="transcript_id",Parent="gene_id"))

  # read the gff3/gtf file
  dat <- readGff3(obj,quiet = !verbose)

  # Check that all required features are there
  min.input <- switch(input,
            "gff3"=c("gene","exon"),
            "gtf"=c("exon"))

  if(!all(min.input %in% levels(dat$type))){
    stop("Your annotation file does not contain all the necessary features. It should contain 'gene' (for gff3 only), 'exon' and at least one of the 'features' argument")
  }

  # If gtf, reformat the attributes and drop the double quotes
  if(input=="gtf"){
    dat$gffAttributes <- gsub("\"","",.convertGffToGtfAttributes(dat$gffAttributes))
  }

  ## get the gene <-> mRNA/transcript map
  # This is mRNA IDs and their parents (genes)
  sel <- dat$type %in% features

  # That step would not necessary for gtf, but it is easier to implement in a
  # similar way for both format
  idMap <- data.frame(type = dat[sel]$type,
                      getGffAttribute(dat[sel],relation$ID),
                      getGffAttribute(dat[sel],relation$Parent))

  ## extract the exons and group by gene ID
  exon.sel <- dat$type == "exon"

  # If we have a minimal gtf (no gene); get the idMap
  if( input=="gtf" & nrow(idMap) ==0 ){
      idMap <- data.frame(type = dat[exon.sel]$type,
                          getGffAttribute(dat[exon.sel],relation$ID),
                          getGffAttribute(dat[exon.sel],relation$Parent))
  }

  # fail if we get nothing - this should not happen though
  stopifnot(nrow(idMap)>0)

  ## we can drop multiple Parents (i.e. comma separated Parent values as we are
  ## collapsing them anyway)
  mRnaID <- sub(",.*","",getGffAttribute(dat[exon.sel],switch(input,
                                                         "gff3"=relation$Parent,
                                                         "gtf"=relation$ID)))

  ## avoid unwanted features
  rngs <- IRanges(start = dat[exon.sel, 1],
                  end = dat[exon.sel, 2])[mRnaID %in% idMap[,relation$ID]]

  ## create a set of synthetic exons
  rngList <- reduce(
    split(rngs,
          idMap[match(mRnaID[mRnaID %in% idMap[,relation$ID]],
                      idMap[,relation$ID]),relation$Parent]))

  ## export the gene, exon and features as gff3
  ## create the new gff object
  ## select the gene
  gene.sel <- dat$type == "gene"

  ## create the gene gff
  geneID <- getGffAttribute(dat[gene.sel,],
                            switch(input,
                                   "gff3"=relation$ID,
                                   "gtf"=relation$Parent))
  geneGff <- switch(input,
                    "gff3"={dat[gene.sel][geneID %in%
                                              idMap[,relation$Parent]]},
                    "gtf"={
                        if (sum(gene.sel)==0){
                            # Note: this is not robust against gene spanning
                            # both strands - get the first exon
                            geneGff <- dat[exon.sel,][match(
                                unique(idMap[idMap$type=="exon",relation$Parent]),
                                idMap[idMap$type=="exon",relation$Parent]),]
                            # convert the type
                            geneGff$type <- "gene"
                            # and define the coords
                            ranges <- range(rngList)
                            IDs <- getGffAttribute(geneGff,relation$Parent)[,1]
                            pos <- match(IDs,names(ranges))
                            geneGff[,1] <- unlist(start(ranges))[pos]
                            geneGff[,2] <- unlist(end(ranges))[pos]
                            geneGff$gffAttributes <- paste0("ID=",IDs)
                        } else {
                            geneGff <- dat[gene.sel][geneID %in%
                                                    idMap[,relation$Parent]]
                            geneGff$gffAttributes <- sub(relation$Parent,
                                                         "ID",geneGff$gffAttributes)
                        }
                        geneGff
                    })

  # define the gene ID if gtf is exon only
  if(input=="gtf" & sum(gene.sel)==0){
      geneID <- getGffAttribute(geneGff,"ID")
  }

  ## create gffs for each feature
  feats <- lapply(features, function(f) {
    if(verbose){message(sprintf("Processing %s features",f))}
    f.sel <- geneID %in% idMap[,relation$Parent][idMap$type == f]
    fGff <- NULL
    if(sum(f.sel)>0){
      fGff <- dat[sel][f.sel]
      fGff$type <- f
      fGff$gffAttributes <- paste("ID=",
                                  getGffAttribute(fGff,relation$ID),
                                  ".0;Parent=",
                                  getGffAttribute(fGff,relation$ID),
                                  sep="")
    }
    return(fGff)
  })
  featureGff <- Reduce(c,feats[!sapply(feats,is.null)])

  # create the featureGff if we have only exon gtf
  if(input=="gtf" & is.null(featureGff)){
      featureGff <- geneGff
      featureGff$type <- "mRNA"
      featureGff$gffAttributes <- paste("ID=",
                          geneID,
                          ".0;Parent=",
                          geneID,
                          sep=""
      )
  }

  ## create the exon gff
  rngList <- rngList[match(geneID[geneID %in% idMap[,relation$Parent]], names(rngList))]
  exonNumber <- elementLengths(rngList)
  if(sum(gene.sel)>0){
      exonGff <- dat[rep(which(gene.sel)[geneID %in% idMap[,relation$Parent]], exonNumber)]
  } else {
      exonGff <- geneGff[rep(1:length(exonNumber),exonNumber)]
  }
  exonGff[,1] <- unlist(start(rngList))
  exonGff[,2] <- unlist(end(rngList))

  exonID <- sapply(exonNumber, ":", 1)
  sel <- geneGff$strand == "+"
  exonID[sel] <- sapply(exonID[sel], rev)
  ID <- getGffAttribute(exonGff, switch(input,
                                        "gff3"=relation$ID,
                                        "gtf"={
                                            ifelse(sum(gene.sel)==0,
                                                   "ID",
                                            relation$Parent)}))
  exonGff$gffAttributes <- paste0("ID=",
                                  paste(ID,
                                        "exon",
                                        unlist(exonID, use.names=FALSE),
                                        sep="."),
                                  ";Name=",
                                  paste(ID,
                                        "exon",
                                        unlist(exonID, use.names=FALSE),
                                        sep="."),
                                  ";Parent=",
                                  paste(ID,"0",sep = "."))
  exonGff$type <- "exon"

  ## combine
  newgff <- c(geneGff, featureGff, exonGff)

  ## change the source
  newgff$source <- "easyRNASeq"

  ## sort
  newgff <- newgff[order(seqnames(newgff), newgff[, 1],
                         factor(as.character(newgff$type),
                                labels = seq_len(2 + length(features)),
                                levels = c("gene", features, "exon"))), ]

  return(switch(output,
                "Genome_intervals" = newgff,
                "GRanges" = {
                  grng <- as(newgff[newgff$type == "exon"], "GRanges")
                  elt <- elementMetadata(grng)[,colnames(elementMetadata(grng))
                                               %in% c("ID","Parent")]
                  colnames(elt) <- c("exon","transcript")
                  elt$gene <- sub("\\.0","",elt$transcript)
                  elementMetadata(grng) <- elt
                  grng
                }
                ))
})
