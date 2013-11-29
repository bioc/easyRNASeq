## TODO make sure to provide the possibility to load a gene model
## test the code: sapply(dir("/Users/delhomme/Documents/EMBL/Projects/HTFGC/NGS/trunk/src/R/packages/easyRNASeq/R",full.names=TRUE),source)

## TODO we need to unify the exons/features, i.e. make sure they are unique
## the easiest way is probably to report the unique only whenever exons or
## features are used. for transcripts and genes, we need to ensure that this is the case

## TODO check why the lib size are different when calculated by edgeR and by RNAseq

## TODO think of using match.arg for default values
## match.arg(c("auto", "variableStep", "fixedStep"),c("auto", "variableStep", "fixedStep"))
## and to replace the .checkArguments function actually!!


## get the annotation
## TODO check the GenomicFeatures package

##' Compute the coverage from a Short Read Alignment file
##' 
##' Computes the genomic reads' coverage from a
##' read file in bam format or any format supported by \pkg{ShortRead}.
##' 
##' \dots{} for fetchCoverage: Can be used for readAligned method from package
##' \pkg{ShortRead} or for scanBamFlag method from package \pkg{Rsamtools}.
##' 
##' @aliases fetchCoverage
##' @name easyRNASeq coverage methods
##' @rdname easyRNASeq-coverage-methods
##' @param obj An \code{\linkS4class{RNAseq}} object
##' @param chr.map A data.frame describing the mapping of original chromosome
##' names towards wished chromosome names. See details.
##' @param chr.sel A vector of chromosome names to subset the final results.
##' @param filename The full path of the file to use
##' @param filter The filter to be applied when loading the data using the
##' "aln" format
##' @param format The format of the reads, one of "aln","bam". If not "bam",
##' all the types supported by the ShortRead package are supported too.
##' @param gapped Is the bam file provided containing gapped alignments?
##' @param ignoreWarnings set to TRUE (bad idea! they have a good reason to be
##' there) if you do not want warning messages.
##' @param isUnmappedQuery additional argument for scanBamFlag \pkg{Rsamtools}
##' @param type The type of data when using the "aln" format. See the
##' \pkg{ShortRead} package.
##' @param validity.check Shall UCSC chromosome name convention be enforced
##' @param what additional argument for ScanBamParam \pkg{Rsamtools}
##' @param \dots additional arguments. See details
##' @return An \code{\linkS4class{RNAseq}} object. The slot readCoverage
##' contains a SimpleRleList object representing a list of coverage vectors,
##' one per chromosome.
##' @author Nicolas Delhomme
##' @seealso \code{\linkS4class{Rle}}
##' \code{\link[ShortRead:readAligned]{ShortRead:readAligned}}
##' @keywords methods
##' @examples
##' 
##' 	\dontrun{
##' 	library("RnaSeqTutorial")
##' 	library(BSgenome.Dmelanogaster.UCSC.dm3)
##' 
##' 	obj <- new('RNAseq',
##' 		organismName="Dmelanogaster",
##' 		readLength=36L,
##' 		chrSize=as.list(seqlengths(Dmelanogaster))
##' 		)
##' 	
##' 	obj <- fetchCoverage(
##' 			obj,
##' 			format="bam",
##'                         filename=system.file(
##' 				"extdata",
##' 				"ACACTG.bam",
##'                             	package="RnaSeqTutorial")
##' 			)
##' 	}
##' 
setMethod(
          f="fetchCoverage",
          signature="RNAseq",
          definition=function(obj,
            format=c("aln","bam"),
            filename=character(1),
            filter=srFilter(),
            type="SolexaExport",
            chr.sel=c(),
            isUnmappedQuery=FALSE,
            what=c("rname","pos","qwidth"),
            validity.check=TRUE,
            chr.map=data.frame(),
            ignoreWarnings=FALSE,
            gapped=TRUE,
            bp.coverage=FALSE,...){
            
            ## check the filename
            if(!file.exists(filename)){
              stop(paste("Cannot read the file:",filename))
            }

            ## set fileName slot if unset (just the filename)
            if(length(fileName(obj)) == 0){
              fileName(obj) <- basename(filename)
            }
            
            ## check the format
            .checkArguments("fetchCoverage","format",format)
            
            ## check if we have the index with bai
            if(format=="bam" & !file.exists(paste(filename,"bai",sep="."))){
              stop(paste("We are missing the index file: ",filename,".bai",sep=""))
            }

            ## are we looking for gapped alignments?
            if(gapped){
              switch(format,
                     "aln"={
                       if(ignoreWarnings){
                         warning("The 'gapped' flag is ignored with data of the 'aln' format kind.")
                       }
                     },
                     "bam"={
                       format="gapped"
                     })
            }
            
            ## switch to read the file
            ## not sure that the ... works for readAligned.
            aln.info <- switch(format,
                               aln=.extractIRangesList(
                                 readAligned(dirname(filename),
                                             pattern=basename(filename),
                                             filter=filter,type=type,...),
                                 chr.sel),
                               bam={
                                 flag <- eval(parse(text=paste(
                                                      "scanBamFlag(isUnmappedQuery=isUnmappedQuery",
                                                      .getArguments(scanBamFlag,...),")",sep="")))
                                 .extractIRangesList(scanBam(filename,
                                                             index=filename,
                                                             param=ScanBamParam(flag=flag,what=what))[[1]],
                                                     chr.sel)
                               },
                               gapped=.extractIRangesList(
                                   readGAlignments(filename,
                                                   index=filename,
                                                   format="BAM",
                                                   use.names=TRUE
                                                   ),
                                 chr.sel
                                 )
                               )
            
            ## stop if the chr sel removes everything!
            if(length(aln.info$rng.list)==0 | sum(elementLengths(aln.info$rng.list)) == 0){
              stop(paste("No data was retrieved from the file: ",
                         filename,
                         ". Make sure that your file is valid, that your 'chr.sel' (if provided) contains valid values; i.e. values as found in the alignment file, not as returned by 'RNAseq'.",
                         sep=""))
            }
            
            librarySize(obj) <- aln.info$lib.size
            
            ## UCSC chr naming convention validity check
            if(validity.check){
              ## modified in version 1.1.9 (06.03.2012) as it was unwise to check for chr in the names
              ## that's dealt with in the .convertToUSCS function
##              chr.grep <- grep("chr",names(aln.info$rng.list))
##              if(length(chr.grep)== 0 | !all(1:length(names(aln.info$rng.list)) %in% chr.grep)){
                if(organismName(obj) != "custom"){
                  if(!ignoreWarnings){
                    warning("You enforce UCSC chromosome conventions, however the provided alignments are not compliant. Correcting it.")
                  }
                }
                names(aln.info$rng.list) <- .convertToUCSC(names(aln.info$rng.list),organismName(obj),chr.map)
##              }
              }

            ## check if we have a single read length
            rL <- unique(do.call("c",lapply(width(aln.info$rng.list),unique)))
            rL <- rL[rL != 0]
            rL <- sort(rL,decreasing=TRUE)
            
            ## check what the user provided
            ## the double && is to make sure we have
            ## a single value tested even if rL
            ## has more than one element. R test anyway all conditions...
            ## 2012-07-06 Changed to change the readLengh as well when there are multiple
            ## readLengths
            if(length(readLength(obj))==1 && readLength(obj)==integer(1)){
              .catn("Updating the read length information.")
              if(length(rL)>1){
                .catn(switch(format,"gapped" = "The alignments are gapped.",
                             "The reads have been trimmed."))
                .catn("Minimum length of",min(rL),"bp.")
                .catn("Maximum length of",max(rL),"bp.")
              } else {                
                .catn("The reads are of",rL,"bp.")
              }
            } else {
              if(any(rL != readLength(obj))){
                warning(paste("The read length stored in the object (probably provided as argument):",
                              readLength(obj),
                              "\nis not the same as the",ifelse(length(rL)==1,"one:","ones:"),
                              paste(rL,collapse=", "),"determined from the file:",
                              filename,"\nUpdating the readLength slot."))
              }
            }
            readLength(obj) <- as.integer(rL)
            
            ## check for the chromosome size and report any problem
            ## TODO this would not be necessary if chr.size is auto
            if(!all(names(aln.info$rng.list) %in% names(chrSize(obj)))){
              warning(paste("The chromosome(s):",
                            paste(names(aln.info$rng.list)[!names(aln.info$rng.list) %in% names(chrSize(obj))],collapse=", "),
                            "is (are) not present in the provided 'chr.sizes' argument"))
            }
           
            if(any(unlist(start(aln.info$rng.list),use.names=FALSE) > chrSize(obj)[rep(names(aln.info$rng.list),
                                                     elementLengths(start(aln.info$rng.list)))])){
              stop("Some of your read coordinates are bigger than the chromosome sizes you provided. Aborting!")
            }
            
            ## check and correct the names in the width and in the ranges, keep the common selector
            valid.names <- sort(intersect(names(aln.info$rng.list),names(chrSize(obj))))
            if(length(chr.sel)>0){
              chrs <- .convertToUCSC(chr.sel,organismName(obj),chr.map)
              if(!all(chrs %in% valid.names)){
                valid.names <- valid.names[valid.names %in% chrs]
                if(!ignoreWarnings){
                  warn=FALSE
                  if(!all(names(aln.info$rng.list)[names(aln.info$rng.list) %in% chrs] %in% valid.names)){
                    warning("Not all the selected ('chr.sel') chromosome names from your read file(s) (aln or bam) exist in your chromosome size list 'chr.sizes'.")   
                    warn=TRUE
                  }
                  if(!all(names(chrSize(obj))[names(chrSize(obj)) %in% chrs] %in% valid.names)){
                    warning("Not all the selected ('chr.sel') chromosome names from the chromosome size list 'chr.sizes' are present in your read file(s) (aln or bam).")
                    warn=TRUE
                  }
                  if(warn & !ignoreWarnings){
                    warning(paste("The available chromosomes in both your read file(s) (aln or bam) and 'chr.sizes' list were restricted to their common term.\n",
                                  "These are: ",paste(valid.names,collapse=", "),".",sep=""))
                  }
                }
              }
            } else {
              if(!ignoreWarnings){
                warn=FALSE
                if(!all(names(aln.info$rng.list) %in% valid.names)){
                  warning("Not all the chromosome names present in your read file(s) (aln or bam) exist in your chromosome size list 'chr.sizes'.")   
                  warn=TRUE
                }
                if(!all(names(chrSize(obj))%in%valid.names)){
                  warning("Not all the chromosome names in your chromosome size list 'chr.sizes' are present in your read file(s) (aln or bam).")
                  warn=TRUE
                }
                if(warn & !ignoreWarnings){
                  warning(paste("The available chromosomes in both your read file(s) (aln or bam) and 'chr.sizes' list were restricted to their common term.\n",
                                "These are: ",paste(valid.names,collapse=", "),".",sep=""))
                }
              }
            }
            
            ## calc the coverage

            ## define the possibilities
            ## 00 is variable length and read coverage
            ## 01 is variable length and bp coverage
            ## 10 is unique length and read coverage
            ## 11 is unique length and bp coverage
            ## 01 and 11 are the same, we just return the bp coverage

            ## Nico August 6th 2012 v1.3.9
            ## Removed the 1e6 divisions as Herve added the support for numeric Rle as a result
            ## from the coverage vector.
            ## Nico sometime July 2012 v1.3.7
            ## the 1e6 series of division allows us to take into account the read proportions for
            ## read of variable length. One would normally divide by the individual readLength, but
            ## weight can only take integers. As a consequence, using a multiplier / divisor of 1e6
            ## lessen the effect of rounding up
            ## IMPORTANT: note that this might result in an integer overflow in the coverage function
            ## that is not reported!! On a 32bit machine, we should still be able to deal with a 2000+ bp coverage
            ## would anyone sequence that deep??
            readCoverage(obj) <- switch(paste(as.integer(c(length(rL)==1,bp.coverage)),collapse=""),
                                        "00" = coverage(aln.info$rng.list[match(valid.names,names(aln.info$rng.list))],
                                          width=as.list(chrSize(obj)[match(valid.names,names(chrSize(obj)))]),
                                          weight=1/width(aln.info$rng.list)[match(valid.names,names(aln.info$rng.list))]),
                                          ##weight=1e6/width(aln.info$rng.list)[match(valid.names,names(aln.info$rng.list))])/1e6,
                                        "10" = coverage(aln.info$rng.list[match(valid.names,names(aln.info$rng.list))],
                                          width=as.list(chrSize(obj)[match(valid.names,names(chrSize(obj)))]))/readLength(obj),
                                        {coverage(aln.info$rng.list[match(valid.names,names(aln.info$rng.list))],
                                                  width=as.list(chrSize(obj)[match(valid.names,names(chrSize(obj)))]))})
          
            ## Nico August 6th 2012 v1.3.9
            ## commented ou as Herve implemented support for Rle numeric vectors for the coverage
            ## and it would return a warning if any coverage value would be NA (i.e. above the 32/64 bit limit). 
            ## Nico sometime July 2012 v1.3.7
            ## ensure that we're not returning junk
            ## could happen if 1e6 is too much and we
            ## reach the integer limits
            ## if(!bp.coverage){
            ##   obs <- sum(sum(readCoverage(obj)))
            ##   exp <- librarySize(obj)
            ##   if( (obs < exp * 0.9) | (obs > exp * 1.1)){
            ##     stop("The observed number of count differs from the expected number! Something went wrong, please contact the author.")
            ##   }
            ## }

            ## Nico August 9th 2012 v.1.3.10
            ## coverage use to return a named list which it does not anymore
            names(readCoverage(obj)) <- valid.names
            
            ## return obj
            return(obj)
          })

## easy call
##' easyRNASeq method
##' 
##' This function is a wrapper around the more low level functionalities of the
##' package.  Is the easiest way to get a count matrix from a set of read
##' files.  It does the following: \itemize{
##' \item{\code{\link[easyRNASeq:ShortRead-methods]{use ShortRead/Rsamtools
##' methods}} for loading/pre-processing the data.}
##' \item{\code{\link[easyRNASeq:fetchAnnotation]{fetch the annotations}}
##' depending on the provided arguments}
##' \item{\code{\link[easyRNASeq:fetchCoverage]{get the reads coverage}} from
##' the provided file(s)}
##' \item{\code{\link[easyRNASeq:easyRNASeq-summarization-methods]{summarize the
##' reads}} according to the selected summarization features}
##' \item{\code{\link[easyRNASeq:easyRNASeq-correction-methods]{optionally
##' apply}} a data correction (i.e. generating RPKM).}
##' \item{\code{\link[easyRNASeq:edgeR-methods]{use edgeR methods}} for
##' post-processing the data or}
##' \item{\code{\link[easyRNASeq:DESeq-methods]{use
##' DESeq methods}} for post-processing the data (either of them being
##' recommended over RPKM).}  }
##' 
##' \itemize{ \item{\dots{} Additional arguments for different functions:
##' \itemize{
##' \item{For the \pkg{biomaRt} \code{\link[biomaRt:getBM]{getBM}} function}
##' \item{For the \code{\link[easyRNASeq:easyRNASeq-annotation-internal-methods]{readGffGtf}}
##' internal function that takes an optional arguments: annotation.type that
##' default to "exon" (used to select the proper rows of the gff or gtf file)}
##' \item{ For the \code{\link[DESeq:estimateDispersions]{DESeq
##' estimateDispersions}} method}
##' \item{For to the \code{\link[base:list.files]{list.files}}
##' function used to locate the read files.}
##' }}
##' \item{the annotationObject When the
##' \code{annotationMethods} is set to \code{env} or \code{rda}, a properly
##' formatted \code{RangedData} or \code{GRangesList} object need to be
##' provided. Check the paragraph RangedData in the vignette or the examples at
##' the bottom of this page for examples. The data.frame-like structure of
##' these objects is where \code{easyRNASeq} will look for the exon, feature,
##' transcript, or gene identifier. Depending on the count method selected, it
##' is essential that the akin column name is present in the annotationObject.
##' E.g. when counting "features", the annotationObject has to contain a
##' "feature" field.}
##' \item{the chr.map The chr.map argument for the easyRNASeq
##' function only works for an "organismName" of value 'custom' with the
##' "validity.check" parameter set to 'TRUE'.  This data.frame should contain
##' two columns named 'from' and 'to'. The row should represent the chromosome
##' name in your original data and the wished name in the output of the
##' function.}
##' \item{count The count can be summarized by exons, features,
##' genes, islands or transcripts. While exons, genes and transcripts are
##' obvious, "features" describes any features provided by the user, e.g.
##' enhancer loci. These are processed as the exons are. For "islands", it is
##' for an under development function that identifies de-novo expression loci
##' and count the number of reads overlapping them. }
##' \item{chr.sizes If set to "auto", then the format has to be "bam", in which
##' case the chromosome names and size are extracted from the BAM header}
##' }
##' 
##' @aliases easyRNASeq easyRNASeq,character-method
##' @rdname easyRNASeq-easyRNASeq
##' @param annotationFile The location (full path) of the annotation file
##' @param annotationObject A \code{\linkS4class{RangedData}} or
##' \code{\linkS4class{GRangesList}} object containing the annotation.
##' @param annotationMethod The method to fetch the annotation, one of
##' "biomaRt","env","gff","gtf" or "rda". All methods but "biomaRt" and "env"
##' require the annotationFile to be set. The "env" method requires the
##' annotationObject to be set.
##' @param chr.map A data.frame describing the mapping of original chromosome
##' names towards wished chromosome names. See details.
##' @param chr.sel A vector of chromosome names to subset the final results.
##' @param chr.sizes A vector or a list containing the chromosomes' size of the
##' selected organism or simply the string "auto". See details.
##' @param conditions A vector of descriptor, each sample must have a
##' descriptor if you use outputFormat DESeq or edgeR. The size of this list
##' must be equal to the number of sample. In addition the vector should be
##' named with the filename of the corresponding samples.
##' @param count The feature used to summarize the reads. One of
##' 'exons','features','genes','islands' or 'transcripts'. See details.
##' @param filenames The name, not the path, of the files to use
##' @param filesDirectory The directory where the files to be used are located.
##' Defaults to the current directory.
##' @param filter The filter to be applied when loading the data using the
##' "aln" format
##' @param format The format of the reads, one of "aln","bam". If not "bam",
##' all the types supported by the \pkg{ShortRead} package are supported too.
##' As of version 1.3.5, it defaults to bam.
##' @param gapped Is the bam file provided containing gapped alignments?
##' @param ignoreWarnings set to TRUE (bad idea! they have a good reason to be
##' there) if you do not want warning messages.
##' @param min.cov When computing read islands, the minimal coverage to take
##' into account for calling an island
##' @param min.length The minimal size an island should have to be kept
##' @param max.gap When computing read islands, the maximal gap size allowed
##' between two islands to merge them
##' @param nbCore defines how many CPU core to use when computing the
##' geneModels. Use the default parallel library
##' @param normalize A boolean to convert the returned counts in RPKM. Valid
##' when the \code{outputFormat} is left undefined (i.e. when a matrix is
##' returned) and when it is \code{DESeq} or \code{edgeR}. Note that it is not
##' advised to normalize the data prior DESeq or edgeR usage!
##' @param organism A character string describing the organism
##' @param outputFormat By default, easyRNASeq returns a matrix.
##' If one of \code{DESeq},\code{edgeR},\code{RNAseq},
##' \code{\linkS4class{SummarizedExperiment}} is provided then
##' the respective object is returned.
##' @param pattern For easyRNASeq, the pattern of file to look for, e.g. "bam$"
##' @param plot Whether or not to plot assessment graphs.
##' @param readLength The read length in bp
##' @param silent set to TRUE if you do not want messages to be printed out.
##' @param summarization A character defining which method to use when
##' summarizing reads by genes. So far, only "geneModels" is available.
##' @param type The type of data when using the "aln" format. See the ShortRead
##' library.
##' @param validity.check Shall UCSC chromosome name convention be enforced?
##' This is only supported for a set of organisms, see
##' \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq:knownOrganisms}},
##' otherwise the argument 'chr.map' can be used to complement it.
##' @param \dots additional arguments. See details
##' @return Returns a count table (a matrix of m features x n samples). If the
##' \code{outputFormat} option has been set, a corresponding object is returned:
##' a \code{\linkS4class{SummarizedExperiment}}, a
##' \code{\link[DESeq:newCountDataSet]{DESeq:newCountDataset}}, a
##' \code{\link[edgeR:DGEList]{edgeR:DGEList}} or \code{\linkS4class{RNAseq}}.
##' @author Nicolas Delhomme
##' @seealso \code{\linkS4class{RNAseq}}
##' \code{\linkS4class{SummarizedExperiment}}
##' \code{\link[edgeR:DGEList]{edgeR:DGEList}}
##' \code{\link[DESeq:newCountDataSet]{DESeq:newCountDataset}}
##' \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq:knownOrganisms}}
##' \code{\link[ShortRead:readAligned]{ShortRead:readAligned}}
##' @keywords methods
##' @examples
##' 
##' 	\dontrun{
##' 	library("RnaSeqTutorial")
##' 	library(BSgenome.Dmelanogaster.UCSC.dm3)
##' 
##' 	## creating a count table from 4 bam files
##' 	count.table <- easyRNASeq(filesDirectory=
##' 		    			system.file(
##' 					"extdata",
##' 					package="RnaSeqTutorial"),
##' 					pattern="[A,C,T,G]{6}\.bam$",
##' 				format="bam",
##' 				readLength=36L,
##' 				organism="Dmelanogaster",
##' 				chr.sizes=as.list(seqlengths(Dmelanogaster)),
##' 				annotationMethod="rda",
##' 				annotationFile=system.file(
##' 				                            "data",
##' 							    "gAnnot.rda",
##' 							    package="RnaSeqTutorial"),
##' 				count="exons")
##' 
##' 	## an example of a chr.map
##' 	chr.map <- data.frame(from=c("2L","2R","MT"),to=c("chr2L","chr2R","chrMT"))
##' 
##' 	## an example of a RangedData annotation
##' 	gAnnot <- RangedData(
##'                      IRanges(
##'                              start=c(10,30,100),
##'                              end=c(21,53,123)),
##'                           space=c("chr01","chr01","chr02"),
##'                           strand=c("+","+","-"),
##'                           transcript=c("trA1","trA2","trB"),
##'                           gene=c("gA","gA","gB"),
##'                           exon=c("e1","e2","e3"),
##'                           universe = "Hs19"
##'                           )
##' 
##' 	## an example of a GRangesList annotation
##' 	grngs <- as(gAnnot,"GRanges")
##' 	grngsList<-split(grngs,seqnames(grngs))
##' }
##' 
## TODO if the summarization ever get changed, modify the if statement when validating the annotation object for no overlapping features
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
            annotationObject = RangedData(),
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

            ## sanity check
            if(!silent){
              .catn("Checking arguments...")
            }

            ## TODO remove in next version
            ## Check if user give a format
            ## if(length(format)>1){
            ##   stop("You must indicate the format of you source files, by setting argument 'format'")
            ## }

            ## we use a default now.
            format <- match.arg(format)

            ## check the chr.sizes
            if(length(chr.sizes)==1){
              if(chr.sizes=="auto" & format != "bam"){
                stop("As you are not using the 'bam' format, you need to set the 'chr.sizes' option.")
              }
              if(is.character(chr.sizes) & chr.sizes != "auto"){
                stop(paste("The 'chr.sizes' option need only be set for non-bam",
                           "formatted files and needs to be a named integer vector."))
              }
            }
            
            ## test the counts
            if(length(count)!=1){
              if(!ignoreWarnings){
                warning("No count method was provided. Defaulting to 'features'.")
              }
              count <- "features"
            }
            .checkArguments("easyRNASeq","count",count)
 
            ## test the summarization
            if(count == "genes" & length(summarization)>1){
              stop(paste("A 'summarization' method is necessary if you choose the 'genes' count method!"))
            }
            
            if(length(summarization)==1){
              .checkArguments("easyRNASeq","summarization",summarization)
            }

            ## check the annotationMethod            
            if(count != "islands"){
              .checkArguments("easyRNASeq","annotationMethod",annotationMethod)
            }
            
            ## check the organism
            if(organism==character(1)){
              if(annotationMethod=="biomaRt"){
                stop("A valid organism name is necessary for the 'organism' arguments when using the 'biomaRt' annotation method.")
              }              
              if(!ignoreWarnings){
                warning("No organism was provided. No validity check for the UCSC compliance of the chromosome name will be applied.")
              }                           
              validity.check=FALSE
            }
            if(!tolower(organism) %in% c(tolower(knownOrganisms()),"custom") & nrow(chr.map) ==0){
              warning(paste("Your organism has no mapping defined to perform the validity check for the UCSC compliance of the chromosome name.",
                            "Defined organism's mapping can be listed using the 'knownOrganisms' function.",
                            "To benefit from the validity check, you can provide a 'chr.map' to your 'easyRNASeq' function call.",
                            "As you did not do so, 'validity.check' is turned off",sep="\n"))
              validity.check=FALSE
            }
            if(organism=="custom" & nrow(chr.map) ==0){
              stop("You want to use a 'custom' organism, but do not provide a 'chr.map'. Aborting.")
            }
            
            ## check the output formats, default to SummarizedExperiments
            ## TODO use a default here as:  format <- match.arg(format)
            outputFormat <- match.arg(outputFormat)

            ## Check if library are loaded
            ## not needed, libraries are loaded by the package
            ## if(0 == length(grep(paste("^package:", 'edgeR',"$", sep=""), search())) & outputFormat=="edgeR"){
            ##   stop("\nLibrary edgeR need to be loaded to use easyRNASeq with option outputFormat equal to 'edgeR'\n")
            ## }
            ## if(0 == length(grep(paste("^package:", 'DESeq',"$", sep=""), search())) & outputFormat=="DESeq"){
            ##   stop("\nLibrary DESeq need to be loaded to use easyRNASeq with option outputFormat equal to 'DESeq'\n")
            ## }
            
            ## check the files
            if((length(filenames) == 0 & pattern == "") | (length(filenames) > 0 & pattern != "")){
              stop("You need to provide either a list of 'filenames' present in the 'filesDirectory' or a 'pattern' matching them.")
            }

            ## if we have filenames, create the pattern
            if(length(filenames) > 0){
              pattern <- paste(filenames, '$',sep="",collapse="|")
            }
            
            ## get source files from the given directory
            filesList <- .list.files(path=path.expand(filesDirectory),pattern=pattern,...)
            names(filesList) <- basename(filesList)
            
            ## check the list of file
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
            
            ## check if we have index with bai
            if(format=="bam"){
              sel <- file.exists(paste(filesList,"bai",sep="."))
              if(any(!sel)){
                stop(paste("Index files (bai) are required. They are missing for the files:",paste(filesList[!sel],collapse = " and ")))
              }
            }

            ## check the conditions
            if(length(conditions)>0){
              if(is.null(names(conditions)) | length(filesList) != length(conditions) | !all(names(filesList) %in% names(conditions))){
                stop("The 'conditions' should be a named vector, the length of the files to proceed. The names should be the names of the files to proceed.")
              }
            }

            ## sort the file lists according to filenames or conditions
            if(length(filenames)>0){
              filesList <- filesList[match(filenames,names(filesList))]
            } else {
              if(length(conditions)>0){
                filesList <- filesList[match(names(conditions),names(filesList))]
              }
            }
            
            ## create the object and fill the fileName
            obj <- new('RNAseq',organismName=organism,readLength=readLength,fileName=names(filesList))
            
            ## Set chromosome size
            if(length(chr.sizes)==1){
              if(chr.sizes == "auto"){
                
                ## read the headers
                headers <- scanBamHeader(filesList)
                
                ## Two sanity checks
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

                ## check if we got some chr sizes at all
                if(length(chr.sizes)==0){
                  stop(paste("No chromosome sizes could be determined from your",
                             "BAM file(s).Is the BAM header present?\nIf not,",
                             "you can use the 'chr.sizes' argument to provide",
                             "the chromosome sizes information."))
                }
              }
            } else {
              if(! is.integer(chr.sizes)){
                stop("chr.sizes should be a named list of integers. 'Use 'as.integer' to convert from numeric.")
              }
              if(is.null(names(chr.sizes))){
                stop("chr.sizes should be a NAMED list of integers. Use 'names()<-' to set the appropriate chromosome names.")
              }
            }

            ## store them
            chrSize(obj) <- chr.sizes
            
            ## check if the chromosome size are valid
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
            
            ## fetch annotation
            if(!silent){
              .catn("Fetching annotations...")
            }
            
            ## provided as an rda?	
            if(annotationMethod=="rda" | annotationMethod == "env"){
              genomicAnnotation(obj) <- switch(
                                               annotationMethod,
                                               "rda" = {
                                                 if(annotationFile==character(1)){
                                                   stop("The annotationMethod 'rda' requires that you provide an 'annotationFile'.")
                                                 }
                                                 if(!file.exists(annotationFile)){
                                                   stop(paste("The provided annotation file:",annotationFile,"does not exist."))
                                                 }
                                                 l.env<-new.env()
                                                 load(annotationFile,envir=l.env)
                                                 if(class(try(gAnnot <- get("gAnnot",envir=l.env)))=="try-error"){
                                                   stop("The provided annotation file does not contain a gAnnot object.")
                                                 }
                                                 if(class(gAnnot) != "RangedData" & class(gAnnot) != "GRangesList"){
                                                   stop("The provided gAnnot object is not of class 'RangedData' or 'GRangesList'")
                                                 }
                                                 gAnnot
                                               },
                                               "env" ={
                                                 if(class(annotationObject) != "RangedData" & class(annotationObject) != "GRangesList"){
                                                   stop("The provided 'annotationObject' object is not of class 'RangedData' or 'GRangesList'")
                                                 }
                                                 if(length(annotationObject)==0){
                                                   stop("The annotationMethod 'env' requires that you provide an 'annotationObject'.")
                                                 }
                                                 annotationObject
                                               })
              
            } else {
              obj <- fetchAnnotation(obj,annotationMethod=annotationMethod,
                                     filename=annotationFile,
                                     ignoreWarnings=ignoreWarnings,...)
            }

            ## check if the annotation contains the valid fields for the count method
            ## check if the annotation are valid
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

              ## check if any annotation is outside the chrSizes boundaries
              common.names <- intersect(names(ranges(obj)),names(chrSize(obj)))
              if(any(sapply(lapply(ranges(obj)[match(common.names,names(ranges(obj)))],range),end) > chrSize(obj)[match(common.names,names(chrSize(obj)))])){
                stop("Your annotation is not in sync with your alignments! Some annotation lie outside the sequences range reported in your BAM file. You may be using two different genome versions.")
              }
              
              ## check for overlaps
              ## TODO this is a bit fishy as it depends on the order of the summarization argument...
              if(!(count == "genes" & summarization[1] == "geneModels")){
                ovl.number <- sum(sapply(findOverlaps(ranges(obj),ignoreSelf=TRUE,ignoreRedundant=TRUE),function(hits){length(unique(queryHits(hits)))}))
                if(ovl.number > 0 & ! ignoreWarnings){
                  warning(paste("There are",ovl.number,"features/exons defined in your annotation that overlap! This implies that some reads will be counted more than once! Is that really what you want?"))
                }
                if(count == "transcripts"){
                  dup.exon <- sum(sapply(findOverlaps(ranges(obj),ignoreSelf=TRUE,type="equal",ignoreRedundant=TRUE),function(hits){length(unique(queryHits(hits)))}))
                  if(dup.exon > 0 & ! ignoreWarnings){
                    warning(paste("There are",dup.exon,"exons defined in your annotation that overlap! This implies that some reads will be counted several time, i.e. once for every transcript! Is that really what you want?"))
                  }
                }
              }
            }
            
            ## check if the chromosome names are valid
            if(validity.check){
              ## TODO what was that for ???
              ## modified in version 1.1.9 (06.03.2012) as it was unwise to check for chr in the names
              ## that's dealt with in the .convertToUSCS function
##              chr.grep <- grep("chr",names(genomicAnnotation(obj)))
##              if(length(chr.grep)== 0 | !all(1:length(names(genomicAnnotation(obj))) %in% chr.grep)){
                if(annotationMethod!="biomaRt" & organismName(obj) != "custom"){
                  if(!ignoreWarnings){
                    warning("You enforce UCSC chromosome conventions, however the provided annotation is not compliant. Correcting it.")
                  }
                }
                ## TODO do I need to put the chr.sel here to ensure we only adapt those selected chromosomes?
                names(genomicAnnotation(obj)) <- .convertToUCSC(names(genomicAnnotation(obj)),organismName(obj),chr.map)
##              }
            }
            
            ## subset the annotation by chr.sel
            if (length(chr.sel) >0){
              if(!chr.sel %in% names(genomicAnnotation(obj))){
                stop(paste("The chromosome name you have given in the 'chr.sel' argument",
                           "does not match any chromosome in your annotation."))
              }
              genomicAnnotation(obj) <- genomicAnnotation(obj)[space(genomicAnnotation(obj)) %in% chr.sel,]
            }
            
            ## Check if the condition list have the same size as the file list
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

            ## Generate the gene model if required
            if(count == 'genes'){
              if(summarization == 'geneModels'){
                if(!silent){
                  .catn("Computing gene models...")
                }
                geneModel(obj) <- .geneModelAnnotation(genomicAnnotation(obj),nbCore)

                ## check the gene model
                ovl.number <- sum(sapply(findOverlaps(geneModel(obj),ignoreSelf=TRUE,ignoreRedundant=TRUE),function(hits){length(unique(queryHits(hits)))}))
                if(ovl.number > 0 & ! ignoreWarnings){
                  warning(paste("There are",ovl.number,"synthetic exons as determined from your annotation that overlap! This implies that some reads will be counted more than once! Is that really what you want?"))
                }
              }
            }
            
            ## Do count
            ## Changed from sapply to lapply to make sure that the rownames are conserved!
            if(!silent){
              .catn("Summarizing counts...")
            }

            ## perform the count (in parallel if asked)
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

            ## decomplex the data
            ## counts
            listOfCount <- do.call(cbind,lapply(countData,function(cData){
              cData$counts
            }))

            ## sizes
            librarySize(obj) <- do.call("c",lapply(countData,function(cData){
              cData$size
            }))
            
            ## we shouldn't get back a list
            if(is.list(listOfCount)){
              warning("Something unexpected happened while calculating the coverage and summarizing it. Aborting and returning the current objects. Check the readCounts slot for more details.")             
              return(list(RNAseq=obj,readCounts=listOfCount))
            }
            
            ## we want proper names!
            colnames(listOfCount) <- fileName(obj)
            if(!all(rownames(listOfCount) %in% .getName(obj,count))){
              warning("Something unexpected happened while calculating the coverage and summarizing it. Aborting and returning the current object. Check the readCounts slot for more details.")
              return(list(RNAseq=obj,readCounts=listOfCount))
            }
            
            ## islands or not
            if( count == 'islands'){
              readCounts(obj)<- .extendCountList(readCounts(obj),listOfCount,count)
            } else{
              readCounts(obj)<- switch(
                                       as.character(length(summarization)),
                                       "1"=.extendCountList(readCounts(obj),listOfCount,count,summarization),
                                       .extendCountList(readCounts(obj),listOfCount,count)
                                       )
            }
               
            ## Return object asked by user
            if(!silent){
              .catn("Preparing output")
            }

            ## if necessary normalize
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
                              ## note that we pass count and summarization as argument to the threedots of the function
                              counts <- .normalizationDispatcher(obj,type="RPKM",count=count,summarization=summarization,plot=FALSE,silent=silent)
                            } else {
                              counts <- readCounts(obj,count,summarization,unique=TRUE)
                            }
                            return(counts)
                          },
                          "SummarizedExperiment"={
                            ## TODO think that for exons/features the count might be redundant
                            if(normalize){
                              if(!ignoreWarnings){
                                warning(paste("Since you want a 'SummarizedExperiment' object,",
                                              "the normalization was not applied to the 'readCounts'",
                                              "slots. Use the RPKM methods on your 'SummarizedExperiment'",
                                              "object to do so."))
                              }
                            }

                            ## get the counts
                            counts <- readCounts(obj,count,summarization)
                            
                            ## create the sample annotation
                            ## TODO should we return the range if we have many reads?
                            ## TODO and at the moment the ReadLength will be 0 if it is not set
                            ## as a parameter... We need to return it as part of the parallel
                            ## processing above - will be easier when all the internals rely on
                            ## SummarizedExperiment
                            colData <- DataFrame(FileName=fileName(obj),
                                                 LibSize=librarySize(obj),
                                                 ReadLength=min(readLength(obj)),
                                                 row.names=fileName(obj))
                            if(length(conditions)>0){
                              colData$Condition <- conditions
                            }

                            ## create the "gene" annotation
                            ## TODO this probably need refactoring...
                            rowData <- switch(count,
                                              "genes"= {
                                                switch(summarization,
                                                       "geneModels"= as(geneModel(obj),"GRanges"),
                                                       switch(class(genomicAnnotation(obj)),
                                                              ## FIXME; need to be tested!
                                                              "RangedData"={
                                                                grng <- as(genomicAnnotation(obj),"GRanges")
                                                                sel <- !duplicated(grng$gene)
                                                                mins <- sapply(split(start(grng),grng$gene),min)
                                                                maxs <- sapply(split(end(grng),grng$gene),max)
                                                                grng <- grng[sel,-match("exon",colnames(grng))]
                                                                start(grng) <- mins[match(grng$gene,names(mins))]
                                                                end(grng) <- maxs[match(grng$gene,names(maxs))]
                                                                if(length(chr.sel)>0){ 
                                                                  seqlevels(grng) <- chr.sel
                                                                  seqnames(grng) <- factor(as.character(seqnames(grng)))
                                                                }
                                                                grng
                                                              },
                                                              ## FIXME; need to be adapted as the above!
                                                              "GRangesList"=unlist(genomicAnnotation(obj)),
                                                              genomicAnnotation(obj)))
                                              },
                                              "transcripts"={
                                                switch(class(genomicAnnotation(obj)),
                                                       "RangedData"={
                                                         grng <- as(genomicAnnotation(obj),"GRanges")
                                                         sel <- !duplicated(grng$transcript)
                                                         mins <- sapply(split(start(grng),grng$transcript),min)
                                                         maxs <- sapply(split(end(grng),grng$transcript),max)
                                                         grng <- grng[sel,-match("exon",colnames(grng))]
                                                         start(grng) <- mins[match(grng$transcript,names(mins))]
                                                         end(grng) <- maxs[match(grng$transcript,names(maxs))]
                                                         if(length(chr.sel)>0){ 
                                                           seqlevels(grng) <- chr.sel
                                                           seqnames(grng) <- factor(as.character(seqnames(grng)))
                                                         }
                                                         grng
                                                         },
                                                       ## FIXME; need to be adapted as the above!
                                                       "GRangesList"=unlist(genomicAnnotation(obj)),
                                                       genomicAnnotation(obj))
                                              },
                                              switch(class(genomicAnnotation(obj)),
                                                     "RangedData"=as(genomicAnnotation(obj),"GRanges"),
                                                     "GRangesList"=unlist(genomicAnnotation(obj)),
                                                     genomicAnnotation(obj))
                                              )
                            ## correct the seq lengths if we have any NA (occurs when we use a RangedData)
                            if(any(is.na(seqlengths(rowData)))){
                              common.names <- intersect(names(chrSize(obj)),names(seqlengths(rowData)))
                              ## no need to check that we have a common set, it was done before so we must have one
                              seqlengths(rowData)[match(common.names,
                                                        names(seqlengths(rowData)))] <- chrSize(obj)[match(common.names,
                                                                                                           names(chrSize(obj)))]
                            }
                            
                            ## the assay contains the data
                            sexp <- SummarizedExperiment(
                                                         assays=SimpleList(counts=counts),
                                                         rowData=rowData,
                                                         colData=colData)

                            ## add the rownames
                            rownames(sexp) <- rownames(counts)

                            ## issue a warning
                            if(any(duplicated(rownames(sexp)))){
                              warning(paste("As you are counting by ",count,", your assay contains redundant entries.",sep=""))
                            }
                            
                            ## report
                            return(sexp)                            
                          }
                          ))
          })

##' count method
##' 
##' This function is to supersed the easyRNASeq function in order to
##' consolidate the option parameters as well as the option output.
##' Ideally, the only output would be a SummarizedExperiment.
##' 
##' @aliases count count,character-method
##' @rdname easyRNASeq-count
##' @param filesDirectory The directory where the files to be used are located.
##' @param outputFormat By default, easyRNASeq returns a
##' \code{\linkS4class{SummarizedExperiment}}. If one
##' of \code{DESeq},\code{edgeR},\code{RNAseq}, \code{matrix} is provided then
##' the respective object is returned. Ideally, this option should get deprecated
##' and only a SummarizedExperiment returned.
##' @param \dots currently additional arguments to the easyRNASeq function.
##' @return Returns a  \code{\linkS4class{SummarizedExperiment}}. If the
##' \code{outputFormat} option has been set, a corresponding object is returned:
##' a count table (a matrix of m features x n samples), a
##' \code{\link[DESeq:newCountDataSet]{DESeq:newCountDataset}}, a
##' \code{\link[edgeR:DGEList]{edgeR:DGEList}} or \code{\linkS4class{RNAseq}}.
##' @author Nicolas Delhomme
##' @seealso \code{\linkS4class{RNAseq}}
##' \code{\linkS4class{SummarizedExperiment}}
##' \code{\link[edgeR:DGEList]{edgeR:DGEList}}
##' \code{\link[DESeq:newCountDataSet]{DESeq:newCountDataset}}
##' \code{\link[easyRNASeq:easyRNASeq-annotation-methods]{easyRNASeq:knownOrganisms}}
##' \code{\link[ShortRead:readAligned]{ShortRead:readAligned}}
##' @keywords methods
##' @examples
##' 
##' 	\dontrun{
##' 	library("RnaSeqTutorial")
##' 	library(BSgenome.Dmelanogaster.UCSC.dm3)
##' 
##' 	## creating a count table from 4 bam files
##' 	sumExp <- count(filesDirectory=system.file(
##'                            "extdata",
##'                            package="RnaSeqTutorial"),
##'                          pattern="[A,C,T,G]{6}\\.bam$",
##'                          readLength=30L,
##'                          organism="Dmelanogaster",
##'                          chr.sizes=seqlengths(Dmelanogaster),
##'                          annotationMethod="rda",
##'                          annotationFile=system.file(
##'                            "data",
##'                            "gAnnot.rda",
##'                            package="RnaSeqTutorial"),
##'                          count="exons"
##'                          )
##'     ## the counts
##'     assays(sumExp)
##'     ## the sample info
##'     colData(sumExp)
##'     ## the 'features' info
##'     rowData(sumExp)
##' }
##' 
setMethod(f="count",
          signature="character",
          definition=function(
            filesDirectory=getwd(),
            outputFormat="SummarizedExperiment",
            ...
            ){
            warning(paste("This function, meant to supersed the 'easyRNASeq' one is under active development.",
                          "Especially the parameters - being consolidated- will be affected.",
                          "Use the 'easyRNASeq' function in a production environment.",sep="\n"))
            easyRNASeq(filesDirectory=filesDirectory,outputFormat=outputFormat,...)
          })
