##' Compute the coverage from a Short Read Alignment file
##' 
##' Computes the genomic reads' coverage from a
##' read file in bam format or any format supported by \pkg{ShortRead}.
##' 
##' \dots{} for fetchCoverage: Can be used for readAligned method from package
##' \pkg{ShortRead}. The use of the dots for the scanBamFlag method from package
##' \pkg{Rsamtools} has been deprecated, as were the 'what' and 'isUnmappedQuery'
##' argument to the function 
##' 
##' @aliases fetchCoverage-deprecated
##' @name easyRNASeq coverage methods
##' @rdname easyRNASeq-coverage-methods
##' @param obj An \code{\linkS4class{RNAseq}} object
##' @param bp.coverage a boolean that default to FALSE to decide whether
##' coverage is to be calculated and stored by bp
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
##' @param paired Is the bam file containing PE reads?
##' @param stranded Is the bam file from a strand specific protocol?
##' @param type The type of data when using the "aln" format. See the
##' \pkg{ShortRead} package.
##' @param validity.check Shall UCSC chromosome name convention be enforced
##' @param ... additional arguments. See details
##' @return An \code{\linkS4class{RNAseq}} object. The slot readCoverage
##' contains a SimpleRleList object representing a list of coverage vectors,
##' one per chromosome.
##' @author Nicolas Delhomme
##' @seealso \code{\linkS4class{Rle}}
##' \code{\link[ShortRead:readAligned]{ShortRead:readAligned}}
##' @keywords methods
##' @examples
##' 
##'   \dontrun{
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
                      validity.check=TRUE,
                      chr.map=data.frame(),
                      ignoreWarnings=FALSE,
                      gapped=TRUE,
                      paired=FALSE,
                      stranded=FALSE,
                      bp.coverage=FALSE,...){
    
    ## warn for deprecation
    if(.getArguments(scanBamFlag,...)!=""){
      .Deprecated(new="fetchCoverage(...)",old="fetchCoverage(...,isUnmappedQuery=FALSE,what=c('rname','pos','qwidth'),...)")
      warning("Passing scanBamFlag arguments has been deprecated. They will be ignored.")
    }
    
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
    
    ## convert a filename to a BamFile
    if(format=="bam" & is.character(filename)){
      filename <- getBamFileList(filename)
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
                                     filter=filter,type=type,withId=TRUE,...),
                         chr.sel),
                       bam={                                 
                         ## flag <- eval(parse(text=paste(
                         ##                      "scanBamFlag(isUnmappedQuery=isUnmappedQuery",
                         ##                      .getArguments(scanBamFlag,...),")",sep="")))
                         .extractIRangesList(
                           .stream(filename),
                           ##   scanBam(filename,
                           ##   index=filename,
                           ##   param=ScanBamParam(flag=flag,what=what))[[1]],
                           chr.sel)
                       },
                       gapped=.extractIRangesList(
                         .stream(filename),
                         chr.sel
                       )
    )
    
    ## stop if the chr sel removes everything!
    if(length(aln.info$rng.list)==0 | sum(elementNROWS(aln.info$rng.list)) == 0){
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
                      path(filename),"\nUpdating the readLength slot."))
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
                                                                               elementNROWS(start(aln.info$rng.list)))])){
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
