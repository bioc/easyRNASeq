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

##' Fetch genic annotation from a gff/gtf file or using biomaRt
##' 
##' The annotation can be retrieved in two ways \itemize{
##' \item{biomaRt}{Use biomaRt and Ensembl to get organism specific annotation.}
##' \item{gff/gtf}{Use a gff or gtf local annotation file.}}
##' When using \pkg{biomaRt}, it is
##' important that the \code{organismName} slot of the
##' \code{\linkS4class{RNAseq}} object is set the prefix of one of the value
##' available using the \pkg{biomaRt}
##' \code{\link[biomaRt:listDatasets]{listDatasets}} function, e.g.
##' "Dmelanogaster".  When reading from a gff/gtf file, a version 3 formatted
##' gff (gtf are modified gff3 from Ensembl) is expected. The function
##' \pkg{genomeIntervals} \code{\link[genomeIntervals:readGff3]{readGff3}} is
##' used to read the data in.
##' 
##' \dots{} are for additional arguments, passed to the \pkg{biomaRt}
##' \code{\link[biomaRt:getBM]{getBM}} function or to the
##' \code{\link[easyRNASeq:easyRNASeq-annotation-internal-methods]{readGffGtf}}
##' internal function that takes an optional arguments: annotation.type that
##' default to "exon". This is used to select the proper rows of the gff or gtf
##' file.
##' 
##' @aliases fetchAnnotation
##' @name easyRNASeq annotation methods
##' @rdname easyRNASeq-annotation-methods
##' @param obj An object of class \code{RNAseq}
##' @param method one of biomaRt, gff, gtf
##' @param filename If the method is gff or gtf, the actual gtf, gff filename
##' @param ignoreWarnings set to TRUE (bad idea! they have a good reason to be
##' there) if you do not want warning messages.
##' @param \dots See details
##' @return A \code{\linkS4class{RangedData}} containing the fetched
##' annotations.
##' @author Nicolas Delhomme
##' @keywords connection data methods
##' @examples
##' 
##' 	\dontrun{
##' 	library("RnaSeqTutorial")
##' 	obj <- new('RNAseq',
##' 		organismName="Dmelanogaster",
##' 		readLength=36L,
##' 		chrSize=as.list(seqlengths(Dmelanogaster))
##' 		)
##' 
##' 	obj <- fetchAnnotation(obj,
##' 				method="gff",
##'                                 filename=system.file(
##' 						"extdata",
##' 						"annot.gff",
##' 						package="RnaSeqTutorial"))
##' 	}
##' 
setMethod(
          f="fetchAnnotation",
          signature="RNAseq",
          definition=function(obj,
            method=c("biomaRt","gff","gtf"),
            filename=character(1),
            ignoreWarnings=FALSE,...){
            
            ## get the methods
            methods <- eval(formals("fetchAnnotation")$method)
            
            ## check the provided one
            if(!method %in% methods){
              stop(paste(
                         "The given method:",
                         method,
                         "is not part of the supported methods:",
                         paste(methods,collapse=", ")))
            }
            
            ## switch depending on the method
            exon.range <- switch(EXPR=method,
                                 "biomaRt"={.getBmRange(organismName(obj),ignoreWarnings=ignoreWarnings,...)},
                                 "gff"={.getGffRange(organismName(obj),filename=filename,ignoreWarnings=ignoreWarnings,...)},
                                 "gtf"={.getGtfRange(organismName(obj),filename=filename,ignoreWarnings=ignoreWarnings,...)}
                                 )
            
            ## update the obj
            genomicAnnotation(obj)<-exon.range
            
            ## return
            return(obj)
          })


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
            gapped=TRUE,...){
            
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
            aln.ranges <- switch(format,
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
                                   readGappedAlignments(filename,
                                                        index=filename,
                                                        format="BAM"
                                                        ),
                                   chr.sel
                                   )
                                 )
            
            ## stop if the chr sel removes everything!
            if(length(aln.ranges)==0){
              stop(paste("No data was retrieved from the file: ",
                         filename,
                         ". Make sure that your file is valid, that your 'chr.sel' (if provided) contains valid values; i.e. values as found in the file, not as returned by 'RNAseq'.",
                         sep=""))
            } else {
              librarySize(obj) <- sum(as.numeric(sapply(aln.ranges,length)))
            }
            
            ## UCSC chr naming convention validity check
            if(validity.check){
              ## modified in version 1.1.9 (06.03.2012) as it was unwise to check for chr in the names
              ## that's dealt with in the .convertToUSCS function
##              chr.grep <- grep("chr",names(aln.ranges))
##              if(length(chr.grep)== 0 | !all(1:length(names(aln.ranges)) %in% chr.grep)){
                if(organismName(obj) != "custom"){
                  if(!ignoreWarnings){
                    warning("You enforce UCSC chromosome conventions, however the provided alignments are not compliant. Correcting it.")
                  }
                }
                names(aln.ranges) <- .convertToUCSC(names(aln.ranges),organismName(obj),chr.map)
##              }
                
                ## ensure that we have the right readLength and only one length
                rL <- unique(sapply(aln.ranges,function(rng){ifelse(length(rng)>0,unique(width(rng)),0)}))
                rL <- rL[rL != 0]
                if(length(rL) > 1 ){
                  stop(paste("The file", filename, "contains reads of different sizes:",paste(rL,collapse=", "),". We cannot deal with such data at the moment. Please contact the authors to add this functionality." ))
                }
                if(rL != readLength(obj)){
                  warning(paste("The read length stored in the object (probably provided as argument):",
                                readLength(obj),
                                "\nis not the same as the one:",rL,"determined from the file:",
                                filename,"\nUpdating it."))
                  readLength(obj) <- as.integer(rL)
                }
            }

            ## check for the chromosome size and report any problem
            tmp <- sapply(names(aln.ranges),function(chr){
              if(!chr %in% names(chrSize(obj))){
                warning(paste("The chromosome:", chr, "is not present in the provided 'chr.sizes' argument"))
                return(0)
              }
              sum(any(start(aln.ranges[[chr]]) > chrSize(obj)[match(chr,names(chrSize(obj)))]))
            })
            if(any(tmp>0)){
              stop("Some of your read coordinates are bigger than the chromosome sizes you provided. Aborting!")
            }

            ## check and correct the names in the width and in the ranges, keep the common selector
            valid.names <- sort(intersect(names(aln.ranges),names(chrSize(obj))))
            if(length(chr.sel)>0){
              chrs <- .convertToUCSC(chr.sel,organismName(obj),chr.map)
              if(!all(chrs %in% valid.names)){
                valid.names <- valid.names[valid.names %in% chrs]
                if(!ignoreWarnings){
                  warn=FALSE
                  if(!all(names(aln.ranges)[names(aln.ranges) %in% chrs] %in% valid.names)){
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
                if(!all(names(aln.ranges) %in% valid.names)){
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
            readCoverage(obj) <- coverage(aln.ranges[match(valid.names,names(aln.ranges))],width=chrSize(obj)[match(valid.names,names(chrSize(obj)))])
            
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
##' \itemize{ \item{\dots{} Additional arguments, passed
##' to the \pkg{biomaRt} \code{\link[biomaRt:getBM]{getBM}} function or to the
##' \code{\link[easyRNASeq:easyRNASeq-annotation-internal-methods]{readGffGtf}}
##' internal function that takes an optional arguments: annotation.type that
##' default to "exon" (used to select the proper rows of the gff or gtf file)
##' or to the \code{\link[DESeq:estimateDispersions]{DESeq
##' estimateDispersions}} method.}
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
##' selected organism
##' @param conditions A vector of descriptor, each sample must have a
##' descriptor if you use outputFormat DESeq or edgeR. The size of this list
##' must be equal to the number of sample. In addition the vector should be
##' named with the filename of the corresponding samples.
##' @param count The feature used to summarize the reads. One of
##' 'exons','features','genes','islands' or 'transcripts'. See details.
##' @param filenames The name, not the path, of the files to use
##' @param filesDirectory The directory where the files to be used are located
##' @param filter The filter to be applied when loading the data using the
##' "aln" format
##' @param format The format of the reads, one of "aln","bam". If not "bam",
##' all the types supported by the \pkg{ShortRead} package are supported too.
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
##' @param outputFormat By default, easyRNASeq returns a count matrix. If one
##' of \code{DESeq},\code{edgeR},\code{RNAseq} is provided then the respective
##' object will be returned.
##' @param pattern For easyRNASeq, the pattern of file to look for, e.g. "bam$"
##' @param plot Whether or not to plot assessment graphs.
##' @param readLength The read length in bp
##' @param silent set to TRUE if you do not want messages to be printed out.
##' @param summarization A character defining which method to use when
##' summarizing reads by genes. So far, only "geneModels" is available.
##' @param type The type of data when using the "aln" format. See the ShortRead
##' library.
##' @param validity.check Shall UCSC chromosome name convention be enforced
##' @param \dots additional arguments. See details
##' @return Returns a count table (a matrix of m features x n samples) unless
##' the \code{outputFormat} option has been set, in which case an object of
##' type \code{\link[DESeq:newCountDataSet]{DESeq:newCountDataset}} or
##' \code{\link[edgeR:DGEList]{edgeR:DGEList}} or \code{\linkS4class{RNAseq}}
##' is returned
##' @author Nicolas Delhomme
##' @seealso \code{\linkS4class{RNAseq}}
##' \code{\link[edgeR:DGEList]{edgeR:DGEList}}
##' \code{\link[DESeq:newCountDataSet]{DESeq:newCountDataset}}
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
            filesDirectory=character(1),
            organism=character(1),
            chr.sizes=c(),
            readLength=integer(1),
            annotationMethod=c("biomaRt","env","gff","gtf","rda"),
            annotationFile=character(1),
            annotationObject = RangedData(),
            format=c("aln","bam"),
            gapped=FALSE,
            count=c('exons','features','genes','islands','transcripts'),
            outputFormat=c("DESeq","edgeR","matrix","RNAseq"),
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
            
            ## Check if user give a format
            if(length(format)>1){
              stop("You must indicate the format of you source files, by setting argument 'format'")
            }
            .checkArguments("easyRNASeq","format",format)
            
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

            ## check the output formats, default to matrix
            if(length(outputFormat)==4){
              outputFormat='matrix'
            }
            .checkArguments("easyRNASeq","outputFormat",outputFormat)

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
            filesList <- list.files(path.expand(filesDirectory),pattern=pattern,full.name=TRUE)            
            names(filesList) <- basename(filesList)
            
            ## check the list of file
            if(length(filesList) == 0 ){
              stop(
                   paste(
                         "No file to work with, you should check your pattern: '",
                         pattern,
                         "' or your directory:",
                         filesDirectory,
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
                filesList <- filesList[,match(names(conditions),names(filesList))]
              }
            }
            
            ## create the object and fill the fileName
            obj <- new('RNAseq',organismName=organism,readLength=readLength,fileName=names(filesList))
            
            ## Set chromosome size
            if(!is.list(chr.sizes)){
              chr.sizes <- as.list(chr.sizes)
            }
            chrSize(obj) <- chr.sizes
            
            ## check if the chromosome size are valid
            if(validity.check){
##              chr.grep <- grep("chr",names(chrSize(obj)))
##              if(length(chr.grep)== 0 | !all(1:length(names(chrSize(obj))) %in% chr.grep)){
                if(organismName(obj) != "custom"){
                  if(!ignoreWarnings){
                    warning("You enforce UCSC chromosome conventions, however the provided chromosome size list is not compliant. Correcting it.")
                  }
                }
                names(chrSize(obj)) <- .convertToUCSC(names(chrSize(obj)),organismName(obj),chr.map)
##              }
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
              obj <- fetchAnnotation(obj,method=annotationMethod,
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
            countData <- lapply(filesList,function(file,obj=obj,
                                                   format=format,
                                                   filter=filter,
                                                   count=count,
                                                   type=type,
                                                   chr.sel=chr.sel,
                                                   validity.check=validity.check,
                                                   summarization=summarization,
                                                   max.gap=max.gap,
                                                   min.cov=min.cov,
                                                   min.lengh=min.length,
                                                   plot=plot,gapped=gapped,...){
              if(!silent){
                .catn(paste("Processing",basename(file)))
              }
              ## Fetch coverage
              obj <- fetchCoverage(obj,format=format,
                                   filename=file,
                                   filter=filter,type=type,
                                   chr.sel=chr.sel,
                                   validity.check=validity.check,
                                   chr.map=chr.map,
                                   gapped=gapped,...)
              
              ## Do count
              obj <- switch(count,
                            "exons"=exonCounts(obj),
                            "features"=featureCounts(obj),
                            ## no need for the nbCore here, the gene model was already done
                            "genes"=geneCounts(obj,summarization),
                            "transcripts"=transcriptCounts(obj),
                            "islands"=islandCounts(obj,max.gap=max.gap,min.cov=min.cov,min.length=min.length,plot=plot)
                            )

              return(list(counts=readCounts(obj,count,summarization),size=librarySize(obj)))
            },obj,format,filter,count,type,chr.sel,validity.check,summarization,max.gap,min.cov,min.length,plot,gapped,...)

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
                              dgeList <- .normalizationDispatcher(dgeList,type="edgeR",silent=silent)
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
                          }
                          ))
          })

