## the threedots are for:
## RPKM: the count and summarization 
".normalizationDispatcher" <- function(obj,type=c("DESeq","edgeR","RPKM","RNAseq"),silent=FALSE,plot=TRUE,...){

  ## report
  if(!silent){
    .catn("Normalizing counts")
  }
  
  ## get the threedots
  add.args <- list(...)
  
  ## dispatch
  return(switch(type,
                "DESeq"={
                  ## estimate the size factors
                  obj <- estimateSizeFactors( obj )
                  
                  ## calculate the dispersion
                  obj <- estimateDispersions(obj,...)

                  ## plot
                  if(plot){

                    ## this has gotten easy
                    ## the fData holds as many columns as condition or one
                    ## if the method is "pooled" or "blind"
                    ## the dispTable contains the same information
                    ## per condition
                    ## the fitInfo function takes a parameter name
                    ## to return what is needed. This is the name
                    ## of the fData column trimmed of "disp_"
                    switch(as.character(ncol(fData(obj))),
                           "1"=plotDispersionEstimates(obj,cond=sub("disp_","",names(fData(obj)))),
                           {
                             sapply(sub("disp_","",names(fData(obj))),function(nam,obj){
                               plotDispersionEstimates(obj,cond=nam)
                             },obj)
                           })
                  }
                  
                  ## if not silent
                  if(!silent){
                    .catn("The counts have now been normalized by 'DESeq'. You can proceed with you differential expression analysis")
                  }
                  return(obj)
                },
                "edgeR"={
                  ## calculate the size factors
                  obj <- calcNormFactors(obj)

                  if(plot){
                    ## plot the normalization factor per sample pairs
                    apply(combn(rownames(obj$samples),2),2,function(co,obj){plotNormalizationFactors(obj,co[1],co[2])},obj)
                    
                    ## plot the MDS
                    plotMDS.DGEList(obj, main = "MDS of all conditions", labels = rownames(obj$samples))
                  }
                  
                  ## calculate the dispersion
                  obj <- estimateCommonDisp(obj)

                  if(plot){
                    ## plot the dispersion estimate
                    obj <- estimateTagwiseDisp(obj, prior.n = 10, trend = "movingave",prop.used = 0.3, grid.length = 500)
                    plotDispersionEstimates(obj)
                    plotMeanVar(obj, show.raw.vars = TRUE, show.tagwise.vars = TRUE,
                                NBline = TRUE, main="Mean-variance (tag variances against tag abundance)")
                    legend("bottomright",col=c("gray60","lightskyblue","darkred","dodgerblue3",1),
                           pch=c("o","o","x",rep(NA,2)),lty=c(rep(NA,3),1,1),lwd=c(rep(NA,3),4,1),
                           ,pt.cex=0.6,c("raw tagwise variances","gene estimated variance",
                              "100 bin averaged raw variance","common dispersion est. var.",
                              "poisson variance"))
                  }
                  return(obj)
                },
                "RPKM"={

                  ## check them and warn for issues
                  if(! all(names(add.args) %in% c("count","summarization"))){
                    warning("You provided '...' arguments that will not be used!")
                  }

                  if(is.null(add.args$summarization) | is.null(add.args$count)){
                    stop("For getting the normalized matrix from an RNAseq object, you need to provide the 'count' and 'summarization' arguments.")
                  }
                  switch(add.args$count,
                         "genes"={
                           RPKM(obj,from=add.args$summarization,unique=TRUE)
                         },
                         RPKM(obj,from=add.args$count,unique=TRUE)
                         )
                },
                "RNAseq"={
                  warning("Since you want an 'RNAseq' object, there is no normalization to be applied.")
                  obj
                },
                stop(paste("Add the new normalization kind! So far only",
                           paste(sub("matrix","RPKM",eval(formals("easyRNASeq")$outputFormat)),collapse=", "),
                           "are supported."))            
                ))
}

