#' simpleRNASeq method
#'
#' This function is a wrapper around the more low level functionalities of the
#' package. It is the simplest way to get a \code{\linkS4class{RangedSummarizedExperiment}}
#' object from a set of bam files. \code{\linkS4class{RangedSummarizedExperiment}} are
#' containers meant to hold any Next-Generation Sequencing experiment results and
#' metadata. The simpleRNASeq method replaces the
#' \code{\link[easyRNASeq:easyRNASeq-easyRNASeq]{easyRNASeq}} function to
#' simplify the usability. It does the following:
#' \itemize{
#' \item use \code{\link[GenomicAlignments:GAlignments-class]{GenomicAlignments}}
#' for reading/pre-processing the BAM files.
#' \item get the \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{annotations}}
#' depending on the selected parameters
#' \item calculate the coverage from the provided file(s)
#' \item \code{\link[easyRNASeq:easyRNASeq-summarization-methods]{summarizes}} the
#' read counts according to the selected summarization
#' \item returns a \code{\linkS4class{RangedSummarizedExperiment}} object.
#' }
#'
#' @aliases simpleRNASeq simpleRNASeq,BamFileList,RnaSeqParam-method
#' @rdname easyRNASeq-simpleRNASeq
#' @param bamFiles a \code{\linkS4class{BamFileList}} object
#' @param nnodes The number of CPU cores to use in parallel
#' @param override Should the provided parameters override the detected ones
#' @param param RnaSeqParam a \code{\linkS4class{RnaSeqParam}} object
#' that describes the RNA-Seq experimental setup.
#' @param verbose a logical to be report progress or not.
#' @return returns a \code{\linkS4class{RangedSummarizedExperiment}} object.
#' @author Nicolas Delhomme
#' @seealso
#' \itemize{
#' \item{For the input:
#' \itemize{
#' \item \code{\linkS4class{AnnotParam}}
#' \item \code{\linkS4class{BamParam}}
#' \item \code{\linkS4class{RnaSeqParam}}
#' }}
#' \item{For the output:
#' \code{\linkS4class{RangedSummarizedExperiment}}
#' }
#' \item{For related functions:
#' \itemize{
#' \item \code{\linkS4class{BamFile}}
#' \item \code{\linkS4class{BamFileList}}
#' \code{\link[easyRNASeq:easyRNASeq-BamFileList]{getBamFileList}}
#' }
#' }}
#' @keywords methods
#' @examples
#'
#' # the data
#' tdir <- tutorialData()
#' annot <- fetchData("Drosophila_melanogaster.BDGP5.77.with-chr.gtf.gz")
#'
#'  # create the BamFileList, get the BAM and BAI index files from the Bioc cache
#'  filenames <- dir(tdir,pattern="[A,T].*\\.bam$",full.names=TRUE)
#'  indexnames <- sapply(paste0(sub(".*_","",basename(filenames)),".bai"),fetchData)
#'  bamFiles <- getBamFileList(filenames,indexnames)
#'
#'   # create the AnnotParam
#'   annotParam <- AnnotParam(annot,type="gtf")
#'
#'   # create the RnaSeqParam
#'   rnaSeqParam <- RnaSeqParam(annotParam=annotParam)
#'
#'   # get a RangedSummarizedExperiment containing the counts table
#'   sexp <- simpleRNASeq(
#'     bamFiles=bamFiles,
#'     param=rnaSeqParam,
#'     verbose=TRUE
#'   )
#'
#'   # get the counts
#'   assays(sexp)$exons
#'
# TODO integrate that in the right position in the file
# \item{groupBy}{One of genes or chromosomes. The default is "genes". It is the
# prefered way to provide annotation and will results in the use of the Bioconductor
# \code{\link[GenomicRanges:summarizeOverlaps]{summarizeOverlaps}} function. Using the
# "chromosomes" value results in the use of specific
# \code{\link[easyRNASeq:easyRNASeq-count-methods]{easyRNASeq count functionalities}}.}
# TODO implement a way to select the chromosomes
setMethod(f="simpleRNASeq",
          signature=c("BamFileList","RnaSeqParam"),
          definition=function(
            bamFiles=BamFileList(),
            param=RnaSeqParam(),
            nnodes=1,
            verbose=TRUE,
            override=FALSE){

            if(verbose){
              message("==========================")
              message("simpleRNASeq version ",packageVersion("easyRNASeq"))
              message("==========================")
              message("Creating a RangedSummarizedExperiment.")
              message("==========================")
              message("Processing the alignments.")
              message("==========================")
              message("Pre-processing ",length(bamFiles)," BAM files.")
            }

            ## =======================
            # validate the BAMFileList
            ## =======================
            if(verbose){
              message("Validating the BAM files.")
            }
            validate(bamFiles)

            ## =======================
            # create the output object
            ## =======================
            sexp <- SummarizedExperiment(
              colData=DataFrame(
                FilePath=path(bamFiles),
                FileName=basename(names(bamFiles)),
                row.names=basename(names(bamFiles))))

            ## =======================
            # process the bams
            ## =======================
            # 1. the sequence info
            seqinfos <- lapply(bamFiles,seqinfo)
            if(length(seqinfos[[1]])==0){
              stop("It would seem some of your BAM files have no sequence information in their headers. Check these.")
            }
            metadata <- list(SeqInfo=seqinfos[[1]])
            if(!all(sapply(seqinfos,identical,metadata$SeqInfo))){
              stop("Your BAM files do not all have the same sequence informations. Check their headers.")
            }
            metadata(sexp) <- metadata

            if(verbose){
              message("Extracted ",length(metadata$SeqInfo)," reference sequences information.")
            }

            # 2. paired or not and variable length or not
            # TODO edit the yieldSize
            colData(sexp) <- cbind(colData(sexp),
                                   do.call(rbind,
                                           parallelize(bamFiles,.streamForParam,
                                                       nnodes,
                                                       yieldSize=10^6,
                                                       verbose=verbose)))
            # a. paired
            if(verbose){
              message("Found ",sum(!colData(sexp)$Paired)," single-end BAM files.")
              message("Found ",sum(colData(sexp)$Paired)," paired-end BAM files.")
            }

            # TODO can we handle both paired and single - end data at the same time?

            # b. width
            if(verbose){
              dev.null <- sapply(1:length(bamFiles),function(i,rL){
                message("Bam file: ",rownames(rL)[i]," has reads of length ", rL[i,])
              },colData(sexp)[,"ReadLength",drop=FALSE])
            }

            # c. strandedness
            # TODO check the Bioc April 2014 newletter for a method to determine
            # strandedness
            colData(sexp)[,"Stranded"] <- sapply(1:nrow(colData(sexp)),
                                                 function(i,df){
              if(is.na(df[i,"Stranded"])){
                warning(paste("Bam file:",rownames(df)[i],
                              "is considered unstranded."))
                warning(paste("Bam file:",rownames(df)[i],df[i,"Note"]))
                FALSE
              } else {
                #warning("At the moment, this function will improperly count Illumina stranded data, reporting reads that align to the opposite strand of the transcript/gene. Please use the param argument to force an unstranded analysis behaviour. Support for \"reverse\" strand specificity will be implemented ASAP")
                if(verbose){
                  message("Bam file: ",rownames(df)[i],
                          " is ", ifelse(df[i,"Stranded"],
                                         "stranded","unstranded"),".")
                  message(df[i,"Note"])
                }
                df[i,"Stranded"]
              }
            },colData(sexp)[,c("Stranded","Note")])

            # reset the note
            colData(sexp)$Note <- NULL

            # TODO WE need to access the presence of multiple mapping reads


            # d. compare with BamParam if provided.
            # 1. paired
            if(length(unique(colData(sexp)$Paired))!=1){
              warning(paste("You have mixed SE and PE data. Every sample will be",
                            "treated according to the identified characteristic,",
                            "i.e. the provided parameter will be ignored when necessary."))
            } else {
              if(colData(sexp)$Paired[1] != paired(param)){
                warning(paste("You provided an incorrect BAM parameter; ",
                              "'paired' should be set to '",
                              colData(sexp)$Paired[1],
                              "'.",sep=""))
              }
            }

            # update the BAM files
            if(!override){
                if(any(colData(sexp)$Paired)){
                   asMates(bamFiles[colData(sexp)$Paired]) <- TRUE
                }
            } else {
                if(stranded(param)){
                    asMates(bamFiles) <- TRUE
                }
            }

            # 2. stranded
            if(length(unique(colData(sexp)$Stranded))!=1){
              warning(paste("You have mixed stranded and unstranded data. Every sample will be",
                            "treated according to the identified characteristic,",
                            "i.e. the provided parameter will be ignored when necessary."))
            } else {
              if(colData(sexp)$Stranded[1] != stranded(param)){
                if(!override){
                  warning(paste("You provided an incorrect BAM parameter; ",
                                "'stranded' should be set to '",
                                colData(sexp)$Paired[1],
                                "'.",sep=""))
                }
              }
            }

            # use the identified parameters or override
            df <- switch(as.character(override),
                         "FALSE" = colData(sexp)[,c("Paired","Stranded")],
                         "TRUE" = DataFrame(
                           Paired=rep(paired(param),length(bamFiles)),
                           Stranded=rep(stranded(param),length(bamFiles)),
                           row.names=basename(bamFiles)))
            if(override){
              warning(paste("You have chosen to override the detected parameters.",
                            "Hope you know what you are doing.",
                            "Contact me if you think the parameter detection failed."))
            }

            # add the protocol parameter
            warning(paste("As of version 2.15.5, easyRNASeq assumes that, if the",
                          "data is strand specific, the sequencing was done",
                          "using a protocol such as the Illumina TruSeq, where",
                          "the reverse strand is quantified - i.e. the",
                          "strandProtocol argument of the BamParam class defaults",
                          "to 'reverse'."))
            df <- cbind(df,StrandProtocol=rep(strandProtocol(param),length(bamFiles)))

            # 3. not a validation but add the read counts
            colData(sexp)$TotalReads <- parallelize(bamFiles,
                                                    .streamForTotalCount,
                                                    nnodes,verbose=verbose)

            if(verbose){
              dev.null <- sapply(1:length(bamFiles),function(i,tR){
                message("Bam file: ",names(tR)[i]," has ",tR[i]," reads.")
              },colData(sexp)$TotalReads)
            }

            ## =======================
            # validate the annotation
            # and get the annotation
            ## =======================
            if(verbose){
              message("==========================")
              message("Processing the annotation")
              message("==========================")
            }

            # TODO in the annotation validation, we should make sure that
            # the right fields are there and are of the right kind;
            # character or character Rle

            # split into a GRangesList based on transcripts or chromosomes
            # check out why is exon a matrix
            grngs <- getAnnotation(annotParam(param),verbose=verbose)

            # TODO we need to adapt the rda, etc description and the AnnotParam
            # to accept a GRanges and not a GRangesList
            sexp <- as(sexp,"RangedSummarizedExperiment")
            rowRanges(sexp) <- switch(precision(param),
                                    "read"={split(grngs,mcols(grngs)[,sub("s$","",countBy(param))])},
                                    "bp"={split(grngs,seqnames(grngs))})

            if(verbose){
              message("==========================")
              message("Sanity checking")
              message("==========================")
            }

            if(length(intersect(seqnames(metadata(sexp)$SeqInfo),
                                seqnames(seqinfo(sexp))))==0){
              stop("There is no common genomic references between your BAM
                   files and the provided annotation. Fix one or the other.")
            }

            # TODO implement sanity check of chromosome vs chromosome
            # TODO implement sanity check of annotations
            # validate the overlap between the BAM header and the annotation
            # how to validate? check that chromosomes in the BAM file
            # overlap those in the annot
            # what to do if not
            # stop if no overlap
            # stop if not all are present and the option "smth smth" is not set

            ## =======================
            # parallelize
            ## =======================
            if(verbose){
              message("==========================")
              message("Creating the count table")
              message("==========================")
              message("Using ",nnodes, ifelse(nnodes==1," CPU core",
                                              " CPU cores in parallel"))
            }

            # TODO, do as for the HistoneChIPseq package
            # to store the data correctly
            countAssay <- SimpleList(do.call(cbind,
                                             parallelize(bamFiles,
                                                         .streamCount,
                                                         nnodes,
                                                         rowRanges(sexp),
                                                         df,
                                                         param,verbose)))
            names(countAssay) <- countBy(param)
            assays(sexp) <- endoapply(countAssay, unname)

            # done
            if(verbose){
              message("==========================")
              message("Returning a")
              message("      RangedSummarizedExperiment")
              message("==========================")
            }
            return(sexp)
          })
