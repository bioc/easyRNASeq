##' Internal methods of AnnotParam objects
##' 
##' These are \code{\linkS4class{AnnotParam}}{AnnotParam} class internal methods:
##' \itemize{
##' \item .validate validate the content of an AnnotParam object
##' }
##' 
##' @aliases .validate
##' @name AnnotParam internal methods
##' @rdname easyRNASeq-internal-AnnotParam-methods
##' @param obj An AnnotParam object
##' @param verbose To print (or not) messages
##' @param \dots additional arguments passed to the retrieval function.
##' At the moment only forwarded to the \pkg{biomaRt} \code{\link[biomaRt:getBM]{getBM}}
##' function.
##' @return
##' \item{.validate}{invisibly return a TRUE logical on success and stops on failure}
##' @author Nicolas Delhomme
##' @keywords internal

".validate" <- function(obj,verbose=TRUE){
  
  ## for the developer
  stopifnot(is(obj,"AnnotParam"))
  
  ## check that the datasource is set
  if(length(datasource(obj))==0){
    stop("No annotation datasource provided.")
  }
  
  ## we don't check that the type is set, we 
  ## just fail if we don't know it
  ## switch per annotation
  return(switch(
    class(obj),
    "AnnotParamCharacter"=
      switch(
        type(obj),
        "rda"={
          if(!file.exists(datasource(obj))){
            stop(paste("The provided annotation file:",datasource(obj),"does not exist."))
          }      
          l.env<-new.env()
          load(datasource(obj),envir=l.env)
          if(class(try(gAnnot <- get("gAnnot",envir=l.env)))=="try-error"){
            stop("The provided annotation file does not contain a gAnnot object.")
          }
          ## TODO SHOULD NOT WE ALLOW FOR A GRanges instead? getAnnotation returns a GRanges.
          if(class(gAnnot) != "RangedData" & class(gAnnot) != "GRangesList"){
            stop("The provided gAnnot object is not of class 'RangedData' or 'GRangesList'")
          }
          if(verbose){
            message("Validated a datasource of type ",type(obj))
          }
          TRUE
        },    
        "gtf"={
          ## check if the file exists
          if(!file.exists(datasource(obj))){
            stop(paste("The provided annotation file:",datasource(obj),"does not exist."))
          } 
          
          ## read some lines
          ## well we suppose a 1000 is enough to check
          some.lines <- scan(datasource(obj),what=as.list(rep("character",9)),nlines=1000,comment.char="#",quiet=!verbose)
          gffAttr <- do.call(rbind,strsplit(some.lines[[9]]," |;"))
          
          ## stop if the attributes we need are not present
          ## we relax on gene_name
          if(!all(GTF.FIELDS[!GTF.FIELDS %in% c("exon_number","gene_name")] %in% gffAttr[1,])){
            stop(paste("Your gtf file: ",datasource(obj)," does not contain all the required fields: ",
                       paste(GTF.FIELDS[!GTF.FIELDS %in% c("exon_number","gene_name")],collapse=", ")
                       ,".",sep=""))
          }
          if(verbose){
            message("Validated a datasource of type ",type(obj))
          }
          TRUE
        },
        "gff3"={
          ## check if the file exists
          if(!file.exists(datasource(obj))){
            stop(paste("The provided annotation file:",datasource(obj),"does not exist."))
          } 
          
          ## check the header
          if(sub("\\D+","",readLines(datasource(obj),1))!=3){
            stop(paste("Your file:",datasource(obj),
                       "does not contain a gff header: '##gff-version 3' as first line. Is that really a gff3 file?"))
          }
          
          ## read some more and check that we got the proper annotation type
          ## NOTE exon is hardcoded. We might want to change it if we ever
          ## change the .readGffGtf annotation.type parameter
          ## NOTE that we read only 1000 lines including comments and hope it's enough to get a
          ## sufficient validation set
          some.lines <- scan(datasource(obj),
                             what=as.list(rep("character",9)),
                             nlines=1000,comment.char="#",
                             quiet=!verbose,sep="\t")
          if(! all(c("mRNA","exon") %in% some.lines[[3]])){
            stop("The provided gff3 contains no annotation of type 'mRNA' and/or 'exon' in the first 1000 lines.")
          }

          ## select the lines we need
          sel <- some.lines[[3]] %in% c("mRNA","exon")          
          
          ## check for the Parent
          if(length(grep("Parent=",some.lines[[9]][sel]))!=sum(sel)){
            stop("The provided gff3 does not contain a 'Parent' attribute for all the annotation of type 'mRNA' and/or 'exon'.")
          }

          ## check for the ID
          if(length(grep("ID=",some.lines[[9]][sel]))!=sum(sel)){
            stop("The provided gff3 does not contain a 'ID' attribute for all the annotation of type 'mRNA' and/or 'exon'.")
          }
          
          if(verbose){
            message("Validated a datasource of type ",type(obj))
          }
          TRUE      
        },
        "biomaRt"={
          ## check if there is connectivity
          ## check the datasource
          if(length(datasource(obj))==0){
            stop(paste("To use the biomaRt functionnalities, we need a datasource name. Set it using the datasource() function."))
          }
          
          dataset<-paste(tolower(datasource(obj)),"gene_ensembl",sep="_")
          if(! dataset %in% listDatasets(useMart(biomart="ensembl"))$dataset){
            stop(paste("The datasource",datasource,"is not supported by the ensembl biomaRt."))
          }
          if(verbose){
            message("Validated a datasource of type ",type(obj))
          }
          TRUE
        },
        stop(paste("The annotation type",type(obj),"is not implemented."))),
  {
    message("No validation performed at that stage")
    TRUE
  }))
}

".extract" <- function(obj,verbose=TRUE,...){
  
  ## for the developer
  stopifnot(is(obj,"AnnotParam"))
  
  ## switch per annotation
  return(switch(class(obj),
         "AnnotParamCharacter"=switch(
           type(obj),
           "rda"={
             if(verbose){
               message("Retrieving annotation from a ",type(obj)," datasource")
             }
             l.env<-new.env()
             load(datasource(obj),envir=l.env)
             get("gAnnot",envir=l.env)
           },    
          {
            if(verbose){
              message("Retrieving annotation from a ",type(obj)," datasource")
            }
            grngs <- switch(type(obj),
                   "biomaRt"={.getBmRange(obj,...)},
                   "gff3"={.getGffRange(obj)},
                   "gtf"={.getGtfRange(obj)})            
          }),
        {
          if(verbose){
            message("Using the provided annotation as such")
          }
          datasource(obj)
         }))
}
