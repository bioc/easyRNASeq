# TODO we need to unify the exons/features, i.e. make sure they are unique
# the easiest way is probably to report the unique only whenever exons or
# features are used. for transcripts and genes, we need to ensure that this is the case

# TODO think of using match.arg for default values
# match.arg(c("auto", "variableStep", "fixedStep"),c("auto", "variableStep", "fixedStep"))
# and to replace the .checkArguments function actually!!

# easy call
#' easyRNASeq method
#'
#' This function is a wrapper around the more low level functionalities of the
#' package.  Is the easiest way to get a count matrix from a set of read
#' files.  It does the following: \itemize{
#' \item{\code{\link[easyRNASeq:ShortRead-methods]{use ShortRead/Rsamtools
#' methods}} for loading/pre-processing the data.}
#' \item{\code{\link[easyRNASeq:easyRNASeq-annotation-methods]{fetch the annotations}}
#' depending on the provided arguments}
#' \item{\code{\link[easyRNASeq:easyRNASeq-coverage-methods]{get the reads coverage}} from
#' the provided file(s)}
#' \item{\code{\link[easyRNASeq:easyRNASeq-summarization-methods]{summarize the
#' reads}} according to the selected summarization features}
#' \item{\code{\link[easyRNASeq:easyRNASeq-correction-methods]{optionally
#' apply}} a data correction (i.e. generating RPKM).}
#' \item{\code{\link[easyRNASeq:edgeR-methods]{use edgeR methods}} for
#' post-processing the data or}
#' \item{\code{\link[easyRNASeq:DESeq-methods]{use
#' DESeq methods}} for post-processing the data (either of them being
#' recommended over RPKM).}  }
#'
#' \itemize{ \item{\dots{} Additional arguments for different functions:
#' \itemize{
#' \item{For the \pkg{biomaRt} \code{\link[biomaRt:getBM]{getBM}} function}
#' \item{For the \code{\link[easyRNASeq:easyRNASeq-annotation-internal-methods]{readGffGtf}}
#' internal function that takes an optional arguments: annotation.type that
#' default to "exon" (used to select the proper rows of the gff or gtf file)}
#' \item{ For the \code{\link[DESeq:estimateDispersions]{DESeq
#' estimateDispersions}} method}
#' \item{For to the \code{\link[base:list.files]{list.files}}
#' function used to locate the read files.}
#' }}
#' \item{the annotationObject When the
#' \code{annotationMethods} is set to \code{env} or \code{rda}, a properly
#' formatted \code{GRangesList} object need to be
#' provided. Check the vignette or the examples at
#' the bottom of this page for examples. The data.frame-like structure of
#' these objects is where \code{easyRNASeq} will look for the exon, feature,
#' transcript, or gene identifier. Depending on the count method selected, it
#' is essential that the akin column name is present in the annotationObject.
#' E.g. when counting "features", the annotationObject has to contain a
#' "feature" field.}
#' \item{the chr.map The chr.map argument for the easyRNASeq
#' function only works for an "organismName" of value 'custom' with the
#' "validity.check" parameter set to 'TRUE'.  This data.frame should contain
#' two columns named 'from' and 'to'. The row should represent the chromosome
#' name in your original data and the wished name in the output of the
#' function.}
#' \item{count The count can be summarized by exons, features,
#' genes, islands or transcripts. While exons, genes and transcripts are
#' obvious, "features" describes any features provided by the user, e.g.
#' enhancer loci. These are processed as the exons are. For "islands", it is
#' for an under development function that identifies de-novo expression loci
#' and count the number of reads overlapping them. }
#' \item{chr.sizes If set to "auto", then the format has to be "bam", in which
#' case the chromosome names and size are extracted from the BAM header}
#' }
#'
#' @aliases easyRNASeq-defunct easyRNASeq,character-method
#' @rdname easyRNASeq-easyRNASeq
#' @param annotationFile The location (full path) of the annotation file
#' @param annotationObject A
#' \code{\linkS4class{GRangesList}} object containing the annotation.
#' @param annotationMethod The method to fetch the annotation, one of
#' "biomaRt","env","gff","gtf" or "rda". All methods but "biomaRt" and "env"
#' require the annotationFile to be set. The "env" method requires the
#' annotationObject to be set.
#' @param chr.map A data.frame describing the mapping of original chromosome
#' names towards wished chromosome names. See details.
#' @param chr.sel A vector of chromosome names to subset the final results.
#' @param chr.sizes A vector or a list containing the chromosomes' size of the
#' selected organism or simply the string "auto". See details.
#' @param conditions A vector of descriptor, each sample must have a
#' descriptor if you use outputFormat DESeq or edgeR. The size of this list
#' must be equal to the number of sample. In addition the vector should be
#' named with the filename of the corresponding samples.
#' @param count The feature used to summarize the reads. One of
#' 'exons','features','genes','islands' or 'transcripts'. See details.
#' @param filenames The name, not the path, of the files to use
#' @param filesDirectory The directory where the files to be used are located.
#' Defaults to the current directory.
#' @param filter The filter to be applied when loading the data using the
#' "aln" format
#' @param format The format of the reads, one of "aln","bam". If not "bam",
#' all the types supported by the \pkg{ShortRead} package are supported too.
#' As of version 1.3.5, it defaults to bam.
#' @param gapped Is the bam file provided containing gapped alignments?
#' @param ignoreWarnings set to TRUE (bad idea! they have a good reason to be
#' there) if you do not want warning messages.
#' @param min.cov When computing read islands, the minimal coverage to take
#' into account for calling an island
#' @param min.length The minimal size an island should have to be kept
#' @param max.gap When computing read islands, the maximal gap size allowed
#' between two islands to merge them
#' @param nbCore defines how many CPU core to use when computing the
#' geneModels. Use the default parallel library
#' @param normalize A boolean to convert the returned counts in RPKM. Valid
#' when the \code{outputFormat} is left undefined (i.e. when a matrix is
#' returned) and when it is \code{DESeq} or \code{edgeR}. Note that it is not
#' advised to normalize the data prior DESeq or edgeR usage!
#' @param organism A character string describing the organism
#' @param outputFormat By default, easyRNASeq returns a matrix.
#' If one of \code{DESeq},\code{edgeR},\code{RNAseq},
#' \code{SummarizedExperiment} is provided then
#' the respective object is returned.
#' @param pattern For easyRNASeq, the pattern of file to look for, e.g. "bam$"
#' @param plot Whether or not to plot assessment graphs.
#' @param readLength The read length in bp
#' @param silent set to TRUE if you do not want messages to be printed out.
#' @param summarization A character defining which method to use when
#' summarizing reads by genes. So far, only "geneModels" is available.
#' @param type The type of data when using the "aln" format. See the ShortRead
#' library.
#' @param validity.check Shall UCSC chromosome name convention be enforced?
#' This is only supported for a set of organisms, which are
#' Dmelanogaster, Hsapiens, Mmusculus and Rnorvegicus;
#' otherwise the argument 'chr.map' can be used to complement it.
#' @param ... additional arguments. See details
#' @return Returns a count table (a matrix of m features x n samples). If the
#' \code{outputFormat} option has been set, a corresponding object is returned:
#' a \code{\linkS4class{RangedSummarizedExperiment}}, a
#' \code{\link[DESeq:newCountDataSet]{DESeq:newCountDataset}}, a
#' \code{\link[edgeR:DGEList]{edgeR:DGEList}} or \code{\linkS4class{RNAseq}}.
#' @author Nicolas Delhomme
#' @seealso \code{\linkS4class{RNAseq}}
#' \code{\linkS4class{RangedSummarizedExperiment}}
#' \code{\link[edgeR:DGEList]{edgeR:DGEList}}
#' \code{\link[DESeq:newCountDataSet]{DESeq:newCountDataset}}
#' \code{\link[ShortRead:readAligned]{ShortRead:readAligned}}
#' @keywords methods
#' @examples
#'  \dontrun{
#' 	library(BSgenome.Dmelanogaster.UCSC.dm3)
#'
#'  # get the example data files - we retrieve a set of example bam files
#'  # from GitHub using curl, as well as their index.
#'  invisible(sapply(c("ACACTG","ACTAGC"),function(bam){
#'      download.file(paste0("https://github.com/UPSCb/UPSCb/raw/",
#'                           "master/tutorial/easyRNASeq/",bam,".bam"),paste0(bam,".bam"))
#'      download.file(paste0("https://github.com/UPSCb/UPSCb/raw/",
#'                           "master/tutorial/easyRNASeq/",bam,".bam.bai"),paste0(bam,".bam.bai"))
#'  }))
#'
#'  # get an example annotation file - we retrieve it from GitHub using curl
#'  invisible(download.file(paste0("https://github.com/UPSCb/UPSCb/raw/",
#'        "master/tutorial/easyRNASeq/gAnnot.rda"),"gAnnot.rda"))
#'
#' 	# creating a count table from 4 bam files
#' 	count.table <- easyRNASeq(filesDirectory=".",
#' 					pattern="[A,C,T,G]{6}\\.bam$",
#' 				format="bam",
#' 				readLength=36L,
#' 				organism="Dmelanogaster",
#' 				chr.sizes=seqlengths(Dmelanogaster),
#' 				annotationMethod="rda",
#' 				annotationFile="gAnnot.rda",
#' 				count="exons")
#'
#' 	# an example of a chr.map
#' 	chr.map <- data.frame(from=c("2L","2R","MT"),to=c("chr2L","chr2R","chrMT"))
#'
#' 	# an example of a GRangesList annotation
#' 	grngs <- GRanges(seqnames=c("chr01","chr01","chr02"),
#'                      ranges=IRanges(
#'                              start=c(10,30,100),
#'                              end=c(21,53,123)),
#'                           strand=c("+","+","-"),
#'                           transcript=c("trA1","trA2","trB"),
#'                           gene=c("gA","gA","gB"),
#'                           exon=c("e1","e2","e3")
#'                           )
#'
#' 	grngsList<-split(grngs,seqnames(grngs))
#' }
setMethod(
    f="easyRNASeq",
    signature="character",
    definition=function(
        filesDirectory=getwd(),
        organism=character(1),
        chr.sizes=c("auto"),
        readLength=integer(1),
        annotationMethod=c("biomaRt","env","gff","gtf","rda"),
        annotationFile=character(1),
        annotationObject = GRangesList(),
        format=c("bam","aln"),
        gapped=FALSE,
        count=c('exons','features','genes','islands','transcripts'),
        outputFormat=c("matrix","SummarizedExperiment","DESeq","edgeR","RNAseq"),
        pattern=character(1),filenames=character(0),nbCore=1,
        filter=srFilter(),type="SolexaExport",
        chr.sel=c(),summarization=c("bestExons","geneModels"),
        normalize=FALSE,max.gap=integer(1),min.cov=1L,
        min.length=integer(1),plot=TRUE,
        conditions=c(),validity.check=TRUE,
        chr.map=data.frame(),
        ignoreWarnings=FALSE,
        silent=FALSE,...){

        # deprecation
        .Defunct("simpleRNASeq")

        # sanity check
        if(!silent){
            .catn("Checking arguments...")
        }

        # TODO remove in next version
        # Check if user give a format
        # if(length(format)>1){
        #   stop("You must indicate the format of you source files, by setting argument 'format'")
        # }

        # we use a default now.
        format <- match.arg(format)

        # check the chr.sizes
        if(length(chr.sizes)==1){
            if(chr.sizes=="auto" & format != "bam"){
                stop("As you are not using the 'bam' format, you need to set the 'chr.sizes' option.")
            }
            if(is.character(chr.sizes) & chr.sizes != "auto"){
                stop(paste("The 'chr.sizes' option need only be set for non-bam",
                           "formatted files and needs to be a named integer vector."))
            }
        }

        # test the counts
        if(length(count)!=1){
            if(!ignoreWarnings){
                warning("No count method was provided. Defaulting to 'features'.")
            }
            count <- "features"
        }
        .checkArguments("easyRNASeq","count",count)

        # test the summarization
        if(count == "genes" & length(summarization)>1){
            stop(paste("A 'summarization' method is necessary if you choose the 'genes' count method!"))
        }

        if(length(summarization)==1){
            .checkArguments("easyRNASeq","summarization",summarization)
            if(summarization=="bestExons"){
                .Deprecated(NULL,msg=paste(
                                     "The bestExons summarization is deprecated as its relevance",
                                     "for RNA-Seq counting is unconvincing."))
            } else {
                .Deprecated(NULL,msg=paste(
                                     "Consider using 'synthetic transcripts' as",
                                     "described in the section 7.1 of the vignette instead of the",
                                     "count=genes,summarization=geneModels deprecated paradigm."))
            }
        }

        # check the annotationMethod
        if(count != "islands"){
            .checkArguments("easyRNASeq","annotationMethod",annotationMethod)
        }

        # check the organism
        if(organism==character(1)){
            if(annotationMethod=="biomaRt"){
                stop("A valid organism name is necessary for the 'organism' arguments when using the 'biomaRt' annotation method.")
            }
            if(!ignoreWarnings){
                warning("No organism was provided. No validity check for the UCSC compliance of the chromosome name will be applied.")
            }
            validity.check=FALSE
        }
        # the knowOrganisms is defunct so we just replace it here - only call
        # TAKEN AWAY as the function is defunct and to prevent an R CMD check NOTE
        #knownOrganisms <- eval(formals(easyRNASeq:::.convertToUCSC)$organism)
        # A Placeholder to have the var defined
        knownOrganisms <- ""
        if(!tolower(organism) %in% c(tolower(knownOrganisms),"custom") & nrow(chr.map) ==0){
            warning(paste("Your organism has no mapping defined to perform the validity check for the UCSC compliance of the chromosome name.",
                          "To benefit from the validity check, you can provide a 'chr.map' to your 'easyRNASeq' function call.",
                          "As you did not do so, 'validity.check' is turned off",sep="\n"))
            validity.check=FALSE
        }
        if(organism=="custom" & nrow(chr.map) ==0){
            stop("You want to use a 'custom' organism, but do not provide a 'chr.map'. Aborting.")
        }

        # check the output formats, default to SummarizedExperiment
        outputFormat <- match.arg(outputFormat)

        # check the files
        if((length(filenames) == 0 & pattern == "") | (length(filenames) > 0 & pattern != "")){
            stop("You need to provide EITHER a list of 'filenames' present in the 'filesDirectory' OR a 'pattern' matching them.")
        }

        # if we have filenames, create the pattern
        if(length(filenames) > 0){
            pattern <- paste(filenames, '$',sep="",collapse="|")
        }

        # get source files from the given directory
        filesList <- .list.files(path=path.expand(filesDirectory),pattern=pattern,...)
        names(filesList) <- basename(filesList)

        # check the list of file
        if(length(filesList) == 0 ){
            stop(
                paste(
                    "No file to work with, you should check your pattern: '",
                    pattern,
                    "' or your directory:",
                    paste(filesDirectory,".",sep=""),
                    "Note that no recursive search is performed.",
                    sep=" "
                    )
                )
        }

        # check if we have index with bai
        # actually create a BamFileList
        if(format=="bam"){
            filesList <- getBamFileList(filesList)
        }

        # check the conditions
        if(length(conditions)>0){
            if(is.null(names(conditions)) | length(filesList) != length(conditions) | !all(names(filesList) %in% names(conditions))){
                stop("The 'conditions' should be a named vector, the length of the files to proceed. The names should be the names of the files to proceed.")
            }
        }

        # sort the file lists according to filenames or conditions
        if(length(filenames)>0){
            filesList <- filesList[match(filenames,names(filesList))]
        } else {
            if(length(conditions)>0){
                filesList <- filesList[match(names(conditions),names(filesList))]
            }
        }

        # create the object and fill the fileName
        obj <- new('RNAseq',organismName=organism,readLength=readLength,fileName=names(filesList))

        # Set chromosome size
        if(length(chr.sizes)==1){
            if(chr.sizes == "auto"){

                # read the headers
                headers <- lapply(filesList,scanBamHeader)

                # Two sanity checks
                if(!all(sapply(headers,
                               function(header,expected){
                                   identical(sort(names(header$targets)),expected)
                               },sort(names(headers[[1]]$targets))))){
                    stop("Not all BAM files use the same chromosome names.")
                }
                chr.sizes <- headers[[1]]$targets
                if(!all(sapply(headers, function(header,chr.sizes){
                    identical(header$targets[match(names(chr.sizes),
                                                   names(header$targets))],chr.sizes)
                },chr.sizes))){
                    stop("The chromosome lengths differ between BAM files.")
                }

                # check if we got some chr sizes at all
                if(length(chr.sizes)==0){
                    stop(paste("No chromosome sizes could be determined from your",
                               "BAM file(s).Is the BAM header present?\nIf not,",
                               "you can use the 'chr.sizes' argument to provide",
                               "the chromosome sizes information."))
                }
            }
        } else {
            if(! is.integer(chr.sizes)){
                stop("chr.sizes should be a named integer vector. 'Use 'as.integer' to convert from numeric.")
            }
            if(is.null(names(chr.sizes))){
                stop("chr.sizes should be a NAMED integer vector. Use 'names()<-' to set the appropriate chromosome names.")
            }
        }

        # store them
        chrSize(obj) <- chr.sizes

        # check if the chromosome size are valid
        if(validity.check){
            if(organismName(obj) != "custom"){
                chr.grep <- grep("chr",names(chrSize(obj)))
                if(length(chr.grep)== 0 | !all(1:length(names(chrSize(obj))) %in% chr.grep)){
                    if(!ignoreWarnings){
                        warning("You enforce UCSC chromosome conventions, however the provided chromosome size list is not compliant. Correcting it.")
                    }
                }
            }
            names(chrSize(obj)) <- .convertToUCSC(names(chrSize(obj)),organismName(obj),chr.map)
        }

        # fetch annotation
        if(!silent){
            .catn("Fetching annotations...")
        }

        # validate and get them
        annotParam <- switch(annotationMethod,
                             "biomaRt"=AnnotParam(
                                 datasource=organism,
                                 type=annotationMethod),
                             "env"=AnnotParam(
                                 datasource=annotationObject),
                             AnnotParam(
                                 datasource=annotationFile,
                                 type=annotationMethod))
        genomicAnnotation(obj)<-getAnnotation(annotParam)

        # check if the chromosome names in the annotation are valid
        if(validity.check){
            if(organismName(obj) != "custom"){
                chr.grep <- grep("chr",names(chrSize(obj)))
                if(length(chr.grep)== 0 | !all(1:length(names(chrSize(obj))) %in% chr.grep)){
                    if(!ignoreWarnings){
                        warning("You enforce UCSC chromosome conventions, however the provided chromosome size list is not compliant. Correcting it.")
                    }
                }
            }
            seqlevels(genomicAnnotation(obj)) <- .convertToUCSC(seqlevels(genomicAnnotation(obj)),organismName(obj),chr.map)
            sel <- match(seqlevels(genomicAnnotation(obj)),names(chrSize(obj)))
            seqlengths(genomicAnnotation(obj))[!is.na(sel)] <- na.omit(chrSize(obj)[sel])
        }

        # create
        if(!any(names(seqlengths(genomicAnnotation(obj))) %in% seqnames(obj))){
          stop("There is no common sequence names between your annotation and your BAM!")
        }

        # check if the annotation contains the valid fields for the count method
        # check if the annotation are valid
        if(count != "islands"){
            if(!(sub("s$","",count)) %in% colnames(genomicAnnotation(obj))){
                stop(
                    "The provided annotation does not contain the expected valid column name: ",
                    sub("s$","",count),
                    " for the '",
                    count,"' method.",
                    sep=""
                    )
            }

            # check if any annotation is outside the chrSizes boundaries
            common.names <- intersect(names(ranges(obj)),names(chrSize(obj)))
            if(any(sapply(lapply(ranges(obj)[match(common.names,names(ranges(obj)))],range),end) > chrSize(obj)[match(common.names,names(chrSize(obj)))])){
                stop(paste("Your annotation is not in sync with your alignments!",
                           "Some annotation lie outside the sequences range reported",
                           "in your BAM file. You may be using two different genome versions."))
            }

            # check for overlaps
            # TODO this is a bit fishy as it depends on the order of the summarization argument...
            if(!(count == "genes" & summarization[1] == "geneModels")){
                ovl <- findOverlaps(ranges(obj),
                                    drop.self=TRUE,
                                    drop.redundant=TRUE)
                ovl.number <- sum(sapply(ovl,
                                         function(hits){length(unique(queryHits(hits)))}))
                if(ovl.number > 0){
                    if(! ignoreWarnings){
                        warning(paste("There are",ovl.number,"features/exons defined",
                                      "in your annotation that overlap!",
                                      "This implies that some reads will be counted",
                                      "more than once! Is that really what you want?"))
                    }
                    genomicAnnotation(obj)$overlap <- as.table(ovl) > 0
                }
                if(count == "transcripts"){
                    ovl <- findOverlaps(ranges(obj),drop.self=TRUE,type="equal",
                                        drop.redundant=TRUE)
                    dup.exon <- sum(sapply(ovl,
                                           function(hits){length(unique(queryHits(hits)))}))
                    if(dup.exon > 0){
                        if( ! ignoreWarnings){
                            warning(paste("There are",dup.exon,"exons defined in your",
                                          "annotation that overlap! This implies that",
                                          "some reads will be counted several time, i.e.",
                                          "once for every transcript! Is that really what you want?"))
                        }
                    }
                }
            }
        }

        # check if the chromosome names are valid
        if(validity.check){
          if (annotationMethod != "biomaRt" & organismName(obj) != "custom") {
            chr.grep <- grep("chr", names(genomicAnnotation(obj)))
            if (length(chr.grep) == 0 | !all(1:length(names(genomicAnnotation(obj))) %in% chr.grep)) {
              if (!ignoreWarnings) {
                warning("You enforce UCSC chromosome conventions, however the provided annotation is not compliant. Correcting it.")
              }
              # TAKEN AWAY as the function is defunct and to prevent an R CMD check NOTE
              #names(genomicAnnotation(obj)) <- easyRNASeq:::.convertToUCSC(names(genomicAnnotation(obj)),
              #                                                             organismName(obj), chr.map)
            }
          }
        }

        # subset the annotation by chr.sel
        if (length(chr.sel) >0){
            if(!chr.sel %in% names(genomicAnnotation(obj))){
                stop(paste("The chromosome name you have given in the 'chr.sel' argument",
                           "does not match any chromosome in your annotation."))
            }
            genomicAnnotation(obj) <- genomicAnnotation(obj)[seqnames(genomicAnnotation(obj)) %in% chr.sel]
            seqlevels(genomicAnnotation(obj)) <- seqlevels(genomicAnnotation(obj))[seqlevels(genomicAnnotation(obj)) %in% chr.sel]
        }

        # Check if the condition list have the same size as the file list
        if(outputFormat=="DESeq"|outputFormat=="edgeR" ){
            if(length(conditions)!=length(filesList)){
                stop(paste(
                    "The number of conditions:",
                    length(conditions),
                    "did not correspond to the number of samples:",
                    length(filesList)
                    )
                     )
            }
        }

        # Generate the gene model if required
        if(count == 'genes'){
            if(summarization == 'geneModels'){
                if(!silent){
                    .catn("Computing gene models...")
                }
                geneModel(obj) <- .geneModelAnnotation(genomicAnnotation(obj),nbCore)

                # check the gene model
                ovl <- findOverlaps(geneModel(obj),drop.self=TRUE,drop.redundant=TRUE)
                ovl.number <- sum(sapply(ovl,function(hits){length(unique(queryHits(hits)))}))
                if(ovl.number > 0){
                    if(! ignoreWarnings){
                        warning(paste("There are",ovl.number,"synthetic exons as",
                                      "determined from your annotation that overlap!",
                                      "This implies that some reads will be counted",
                                      "more than once! Is that really what you want?"))
                    }
                    geneModel(obj)$overlap <- as.table(ovl) > 0
                }
            }
        }

        # Do count
        # Changed from sapply to lapply to make sure that the rownames are conserved!
        if(!silent){
            .catn("Summarizing counts...")
        }

        # perform the count (in parallel if asked)
        countData <- parallelize(obj=filesList,fun=.doCount,nnodes=nbCore,
                                 rnaSeq=obj,format=format,
                                 filter=filter,count=count,
                                 type=type,chr.map=chr.map,
                                 chr.sel=chr.sel,
                                 validity.check=validity.check,
                                 summarization=summarization,
                                 max.gap=max.gap,min.cov=min.cov,
                                 min.length=min.length,
                                 plot=plot,gapped=gapped,
                                 silent=silent,...)

        # decomplex the data
        # counts
        listOfCount <- do.call(cbind,lapply(countData,function(cData){
            cData$counts
        }))

        # sizes
        librarySize(obj) <- do.call("c",lapply(countData,function(cData){
            cData$size
        }))

        # we shouldn't get back a list
        if(is.list(listOfCount)){
            warning("Something unexpected happened while calculating the coverage and summarizing it. Aborting and returning the current objects. Check the readCounts slot for more details.")
            return(list(RNAseq=obj,readCounts=listOfCount))
        }

        # we want proper names!
        colnames(listOfCount) <- fileName(obj)
        if(!all(rownames(listOfCount) %in% .getName(obj,count))){
            warning("Something unexpected happened while calculating the coverage and summarizing it. Aborting and returning the current object. Check the readCounts slot for more details.")
            return(list(RNAseq=obj,readCounts=listOfCount))
        }

        # islands or not
        if( count == 'islands'){
            readCounts(obj)<- .extendCountList(readCounts(obj),listOfCount,count)
        } else{
            readCounts(obj)<- switch(
                as.character(length(summarization)),
                "1"=.extendCountList(readCounts(obj),listOfCount,count,summarization),
                .extendCountList(readCounts(obj),listOfCount,count)
                )
        }

        # Return object asked by user
        if(!silent){
            .catn("Preparing output")
        }

        # if necessary normalize
        return(switch(outputFormat,
                      "DESeq"={
                          cds <- newCountDataSet(countData=readCounts(obj,count,summarization,unique=TRUE),conditions=conditions)
                          if(normalize){
                              cds <- .normalizationDispatcher(cds,type="DESeq",silent=silent,plot=plot,...)
                          }
                          return(cds)
                      },
                      "edgeR"={
                          dgeList <- DGEList(counts=readCounts(obj,count,summarization,unique=TRUE),group=conditions)
                          if(normalize){
                              dgeList <- .normalizationDispatcher(dgeList,type="edgeR",silent=silent,plot=plot,...)
                          }
                          return(dgeList)
                      },
                      "RNAseq"={
                          if(normalize){
                              if(!ignoreWarnings){
                                  warning("Since you want an 'RNAseq' object, the normalization was not applied to the 'readCounts' slots. Use the RPKM methods on your 'RNAseq' object to do so.")
                              }
                          }
                          return(obj)
                      },
                      "matrix"= {
                          if(normalize){
                              # note that we pass count and summarization as argument to the threedots of the function
                              counts <- .normalizationDispatcher(obj,type="RPKM",count=count,summarization=summarization,plot=FALSE,silent=silent)
                          } else {
                              counts <- readCounts(obj,count,summarization,unique=TRUE)
                          }
                          return(counts)
                      },
                      "SummarizedExperiment"={
                          # TODO think that for exons/features the count might be redundant
                          if(normalize){
                              if(!ignoreWarnings){
                                  warning(paste("Since you want a 'RangedSummarizedExperiment' object,",
                                                "the normalization was not applied to the 'readCounts'",
                                                "slots. Use the RPKM methods on your 'RangedSummarizedExperiment'",
                                                "object to do so."))
                              }
                          }

                          # get the counts
                          counts <- readCounts(obj,count,summarization)

                          # create the sample annotation
                          # TODO should we return the range if we have many reads?
                          # TODO and at the moment the ReadLength will be 0 if it is not set
                          # as a parameter... We need to return it as part of the parallel
                          # processing above - will be easier when all the internals rely on
                          # SummarizedExperiment
                          colData <- DataFrame(FileName=fileName(obj),
                                               LibSize=librarySize(obj),
                                               ReadLength=min(readLength(obj)),
                                               row.names=fileName(obj))
                          if(length(conditions)>0){
                              colData$Condition <- conditions
                          }

                          # create the "gene" annotation
                          # TODO this probably need refactoring...
                          rowRanges <- switch(count,
                                            "genes"= {
                                                switch(summarization,
                                                       "geneModels"= as(geneModel(obj),"GRanges"),
                                                       switch(class(genomicAnnotation(obj)),
                                                              # FIXME; need to be tested!
                                                              # "RangedData"={
                                                              #     grng <- as(genomicAnnotation(obj),"GRanges")
                                                              #     sel <- !duplicated(grng$gene)
                                                              #     mins <- sapply(split(start(grng),grng$gene),min)
                                                              #     maxs <- sapply(split(end(grng),grng$gene),max)
                                                              #     grng <- grng[sel,-match("exon",colnames(grng))]
                                                              #     start(grng) <- mins[match(grng$gene,names(mins))]
                                                              #     end(grng) <- maxs[match(grng$gene,names(maxs))]
                                                              #     if(length(chr.sel)>0){
                                                              #         seqlevels(grng) <- chr.sel
                                                              #         seqnames(grng) <- factor(as.character(seqnames(grng)))
                                                              #     }
                                                              #     grng
                                                              # },
                                                              # FIXME; need to be adapted as the above!
                                                              "GRangesList"=unlist(genomicAnnotation(obj)),
                                                              genomicAnnotation(obj)))
                                            },
                                            "transcripts"={
                                                switch(class(genomicAnnotation(obj)),
                                                       # "RangedData"={
                                                       #     grng <- as(genomicAnnotation(obj),"GRanges")
                                                       #     sel <- !duplicated(grng$transcript)
                                                       #     mins <- sapply(split(start(grng),grng$transcript),min)
                                                       #     maxs <- sapply(split(end(grng),grng$transcript),max)
                                                       #     grng <- grng[sel,-match("exon",colnames(grng))]
                                                       #     start(grng) <- mins[match(grng$transcript,names(mins))]
                                                       #     end(grng) <- maxs[match(grng$transcript,names(maxs))]
                                                       #     if(length(chr.sel)>0){
                                                       #         seqlevels(grng) <- chr.sel
                                                       #         seqnames(grng) <- factor(as.character(seqnames(grng)))
                                                       #     }
                                                       #     grng
                                                       # },
                                                       # FIXME; need to be adapted as the above!
                                                       "GRangesList"=unlist(genomicAnnotation(obj)),
                                                       genomicAnnotation(obj))
                                            },
                                            switch(class(genomicAnnotation(obj)),
                                                   "GRangesList"=unlist(genomicAnnotation(obj)),
                                                   genomicAnnotation(obj))
                                            )

                          # REMOVE?: RangedData are deprecated, keep to check if this is still needed
                          # correct the seq lengths if we have any NA (occurs when we use a RangedData)
                          # if(any(is.na(seqlengths(rowRanges)))){
                          #     common.names <- intersect(names(chrSize(obj)),names(seqlengths(rowRanges)))
                          #     # no need to check that we have a common set, it was done before so we must have one
                          #     seqlengths(rowRanges)[match(common.names,
                          #                               names(seqlengths(rowRanges)))] <- chrSize(obj)[match(common.names,
                          #                                                                                  names(chrSize(obj)))]
                          # }

                          # the assay contains the data
                          sexp <- SummarizedExperiment(
                              assays=SimpleList(counts=counts),
                              rowRanges=rowRanges,
                              colData=colData)

                          # add the rownames
                          rownames(sexp) <- rownames(counts)

                          # issue a warning
                          if(any(duplicated(rownames(sexp)))){
                              warning(paste("As you are counting by ",count,", your assay contains redundant entries.",sep=""))
                          }

                          # report
                          return(sexp)
                      }
                      ))
    })

