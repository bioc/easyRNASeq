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

## calculate the coverage
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

              ## emergency stop
              if(length(intersect(names(readCoverage(obj)),names(genomicAnnotation(obj))))==0 & organism == character(1) & validity.check == FALSE){
                stop(paste("Emergency stop.",
                           "The chromosome names in your bam file do not match those in your annotation.",
                           "You might solve that issue by providing a value to the 'organism' parameter and",
                           "making sure that the 'validity.check' is set to 'TRUE'."))
              }
              
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

