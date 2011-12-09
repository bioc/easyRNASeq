## Extend ShortRead

## an additional filter
chastityFilter <- function(.name="Illumina Chastity Filter") 
{
  srFilter(function(x){
    if(any(rownames(varMetadata(alignData(x))) == "filtering")){
      keep<-alignData(x)$filtering=="Y"
    } else {
      warning(paste("The '",.name,"' filter is only valid for Illumina reads.",sep=""))
      keep<-rep(TRUE,length(x))
    }
    return(keep)
  },name=.name)
}

## de-multiplex multiplexed libs
setMethod(
          f="demultiplex",
          signature="AlignedRead",
          definition=function(obj,barcodes=c(),barcodes.qty=12,barcode.length=6, edition.dist=2, type=c("independant","within"),index.only=FALSE){

            ## TODO we only want one type!!
            ## default to independant
            
            ## check the input
            types <- eval(formals("demultiplex")$type)
            if(!type %in% types){
              stop(paste(
                         "The given type:",
                         type,
                         "is not part of the supported types:",
                         paste(types,collapse=", ")))
            }
            
            ## do we have barcodes
            if(length(barcodes)==0){
              barcodes <- switch(type,
                independant=sort(table(alignData(obj)$multiplexIndex),decreasing=TRUE)[1:barcodes.qty],
                within=sort(table(as.character(narrow(sread(obj),start=1,width=barcode.length))),decreasing=TRUE)[1:barcodes.qty])
            } else {
              if(length(unique(nchar(barcodes)))!=1){
                stop("We accept only barcodes having the same length.")
              }
              if(!all(nchar(barcodes) == barcode.length)){
                warning(paste("The provided barcode length was not correct. Set it to:",nchar(barcodes[1])))
                barcode.length=nchar(barcodes[1])
              }
            }
            
            ## get the barcodes, according to a certain size
            all.barcodes <- switch(type,
                               independant=DNAStringSet(as.character(alignData(obj)$multiplexIndex)),
                               within=narrow(sread(obj),start=1,width=barcode.length)
                               )

            ## just for illumina
            if(type=="independant" & nchar(all.barcodes)[1] > barcode.length){
              all.barcodes <- narrow(all.barcodes,start=1,width=barcode.length)
            }
            
            ## calculate the edition distance
            ## dist<-do.call(cbind,srdistance(barcodes,barcodes))
            ## or this to speed up
            ## calc the dist only on the unique barcodes
            ## rather than calling unique, generate all possibilities instead
            ## and map it back
            unique.barcodes <- DNAStringSet(
                                            apply(
                                                  as.matrix(eval(parse(
                                                                       text=paste(
                                                                         "expand.grid(",
                                                                         paste(rep("c('A','C','G','T','N')",barcode.length),collapse=","),
                                                                         ",stringsAsFactors=FALSE)")
                                                                       ))),1,paste,collapse=""))
            
            dist <- do.call(cbind,srdistance(unique.barcodes,barcodes))[match(all.barcodes,unique.barcodes),]
            
            ## create a selector per barcode
            sels <- lapply(barcodes,
                           function(barcode,dist,edition.dist){
                             i.sel <- colnames(dist)==barcode
                             dist[,i.sel]<=edition.dist & eval(parse(text=paste(paste("dist[,",which(!i.sel),"] > edition.dist",collapse=" & "))))
                           },dist,edition.dist)

            ## validate
            if(sum(sapply(sels,sum))>length(all.barcodes)){
              warning(
                      paste(
                            "You have more reads selected: ",sapply(sels,sum),"using the edition distance:",
                            edition.dist,"than you have reads:",length(all.barcodes),"!\nUse the barcodePlot function to visually assess it."))
            }
            ## select barcodes
            names(sels)<-colnames(dist)

            ## either the alns or the sels are returned
            if(index.only){
              return(sels)
            } else {              
              ## original read length
              read.length <- width(obj)[1]
              
              ## return a lists of aln
              alns <- lapply(barcodes,function(barcode,obj,sels){
                reads <- narrow(obj[sels[[barcode]]],start=barcode.length+1,width=read.length-barcode.length)
                bars <- narrow(obj[sels[[barcode]]],start=1,width=barcode.length)
                return(list(reads=reads,bars=bars))
              },obj,sels)                       
              names(alns) <- colnames(dist)
              alns <- list(reads=lapply(alns,function(x){x[[1]]}),barcodes=lapply(alns,function(x){sread(x[[2]])}))
              return(alns)
            }
          })

setMethod(
          f="demultiplex",
          signature="DNAStringSet",
          definition=function(obj,barcodes=c(),barcodes.qty=12,barcode.length=6, edition.dist=2, type=c("independant","within"),index.only=FALSE){

            ## There's only one possible type!!
            ## TODO extract this error to share it with other methods that go for a DNAStringSet
            if(type=="independant"){
              stop(paste("We cannot accept the independant argument for a ",class(obj),", since we miss the barcode sequences",sep=""))
            }
            
            ## check the input
            types <- eval(formals("demultiplex")$type)[-1]
            if(!type %in% types){
              stop(paste(
                         "The given type:",
                         type,
                         "is not part of the supported types:",
                         paste(types,collapse=", ")))
            }
            
            ## do we have barcodes
            if(length(barcodes)==0){
              barcodes <- sort(table(as.character(narrow(obj,start=1,width=barcode.length))),decreasing=TRUE)[1:barcodes.qty]
            }  else {
              if(length(unique(nchar(barcodes)))!=1){
                stop("We accept only barcodes having the same length.")
              }
              if(!all(nchar(barcodes) == barcode.length)){
                warning(paste("The provided barcode length was not correct. Set it to:",nchar(barcodes[1])))
                barcode.length=nchar(barcodes[1])
              }
            }
            
            ## get the barcodes, according to a certain size
            all.barcodes <-narrow(obj,start=1,width=barcode.length)
            
            ## calculate the edition distance
            ## dist<-do.call(cbind,srdistance(barcodes,barcodes))
            ## or this to speed up
            ## calc the dist only on the unique barcodes
            ## rather than calling unique, generate all possibilities instead
            ## and map it back
            unique.barcodes <- DNAStringSet(
                                            apply(
                                                  as.matrix(eval(parse(
                                                                       text=paste(
                                                                         "expand.grid(",
                                                                         paste(rep("c('A','C','G','T','N')",barcode.length),collapse=","),
                                                                         ",stringsAsFactors=FALSE)")
                                                                       ))),1,paste,collapse=""))
            
            dist <- do.call(cbind,srdistance(unique.barcodes,barcodes))[match(all.barcodes,unique.barcodes),]
            
            ## create a selector per barcode
            sels <- lapply(barcodes,
                           function(barcode,dist,edition.dist){
                             i.sel <- colnames(dist)==barcode
                             dist[,i.sel]<=edition.dist & eval(parse(text=paste(paste("dist[,",which(!i.sel),"] > edition.dist",collapse=" & "))))
                           },dist,edition.dist)

            ## validate
            if(sum(sapply(sels,sum))>length(all.barcodes)){
              warning(
                      paste(
                            "You have more reads selected: ",sapply(sels,sum),"using the edition distance:",
                            edition.dist,"than you have reads:",length(all.barcodes),"!\nUse the barcodePlot function to visually assess it."))
            }
            ## select barcodes
            names(sels)<-colnames(dist)

            ## either the alns or the sels are returned
            if(index.only){
              return(sels)
            } else {              
              ## original read length
              read.length <- width(obj)[1]
              
              ## return a lists of aln
              alns <- lapply(barcodes,function(barcode,obj,sels){
                reads <- narrow(obj[sels[[barcode]]],start=barcode.length+1,width=read.length-barcode.length)
                bars <- narrow(obj[sels[[barcode]]],start=1,width=barcode.length)
                return(list(reads=reads,bars=bars))
              },obj,sels)                       
              names(alns) <- colnames(dist)
              alns <- list(reads=lapply(alns,function(x){x[[1]]}),barcodes=lapply(alns,function(x){x[[2]]}))
              return(alns)
            }
          })

setMethod(
          f="demultiplex",
          signature="ShortReadQ",
          definition=function(obj,barcodes=c(),barcodes.qty=12,barcode.length=6, edition.dist=2, type=c("independant","within"),index.only=FALSE){
            return(
                   demultiplex(sread(obj),barcodes=barcodes,barcodes.qty=barcodes.qty,barcode.length=barcode.length, edition.dist=edition.dist, type=type,index.only=index.only)
                   )
          })

setMethod(
          f="barcodePlot",
          signature="AlignedRead",
          definition=function(obj,barcodes=c(),type=c("independant","within"),barcode.length=6,show.barcode=20,...){

            ## TODO check that we got barcodes
            
            ## check the input
            types <- eval(formals("barcodePlot")$type)
            if(!type %in% types){
              stop(paste(
                         "The given type:",
                         type,
                         "is not part of the supported types:",
                         paste(types,collapse=", ")))
            }
            
            ## get the barcodes
            barcodes <- switch(type,
                               independant=DNAStringSet(as.character(alignData(obj)$multiplexIndex)),
                               within=narrow(sread(obj),start=1,width=barcode.length)
                               )

            ## get the frequency
            freqs <- sort(table(factor(as.character(barcodes)))/length(barcodes),decreasing=TRUE)

            ## get the args
            args <- lapply(as.list(match.call())[-1], eval, parent.frame())
            args <- args[!names(args) %in% rev(names(formals("barcodePlot")))[-1]]

            ## defaults
            x.lim=c(0,1)
            y.lim=c(0,show.barcode)
            p.main=""
            
            ## upd defaults
            if(length(args)>=1){
              ## check if xlim was provided
              if("xlim" %in% names(args)){
                x.lim=rev(1-args$xlim)
              }
              
              ## check if ylim was provided
              if("ylim" %in% names(args)){
                y.lim=args$ylim
              }

              ## check if main was provided
              if("main" %in% names(args)){
                p.main=args$main
              }
            }
            
            ## get remaining args
            args <- args[!names(args)%in%c("xlim","ylim","main")]

            ## plot
            if(length(args)>=1){
              eval(parse(text=paste("plot(0,0,bty='n',type='n',axes=FALSE,xlab='',ylab='',xlim=x.lim,ylim=y.lim,",paste(names(args),args,sep="=",collapse=","),")")))
            } else {
              plot(0,0,bty='n',type='n',axes=FALSE,xlab='',ylab='',xlim=x.lim,ylim=y.lim)
            }
            rect(c(1-freqs[1:show.barcode]),c((show.barcode-1):0),rep(1,show.barcode),c(show.barcode:1),col=ifelse(names(freqs)[1:show.barcode] %in% barcodes,"orange","black"))
            axis(2,labels=names(freqs)[1:show.barcode],at=seq(show.barcode-0.5,0.5,-1),las=2,tick=FALSE,line = -2)
            axis(3,labels=rev(seq(0,1,0.5)),at=seq(0,1,0.5))
            mtext("barcodes frequency",side=3,at=0.5,line=3)
            mtext(p.main,side=1,cex=par("cex.main"))
            
            ## return
            invisible(freqs)
          })

## TODO merge that with the above one
## make sure that we have a proper cex for the y axis
setMethod(
          f="barcodePlot",
          signature="DNAStringSet",
          definition=function(obj,barcodes=c(),type=c("independant","within"),barcode.length=6,show.barcode=20,...){

            ## TODO check that we got barcodes
            
            ## check the input
            types <- eval(formals("barcodePlot")$type)
            if(!type %in% types){
              stop(paste(
                         "The given type:",
                         type,
                         "is not part of the supported types:",
                         paste(types,collapse=", ")))
            }
            
            ## get the barcodes
            barcodes <- narrow(obj,start=1,width=barcode.length)

            ## get the frequency
            freqs <- sort(table(factor(as.character(barcodes)))/length(barcodes),decreasing=TRUE)

            ## get the args
            args <- lapply(as.list(match.call())[-1], eval, parent.frame())
            args <- args[!names(args) %in% rev(names(formals("barcodePlot")))[-1]]

            ## defaults
            x.lim=c(0,1)
            y.lim=c(0,show.barcode)
            p.main=""
            
            ## upd defaults
            if(length(args)>=1){
              ## check if xlim was provided
              if("xlim" %in% names(args)){
                x.lim=rev(1-args$xlim)
              }
              
              ## check if ylim was provided
              if("ylim" %in% names(args)){
                y.lim=args$ylim
              }

              ## check if main was provided
              if("main" %in% names(args)){
                p.main=args$main
              }
            }
            
            ## get remaining args
            args <- args[!names(args)%in%c("xlim","ylim","main")]

            ## plot
            if(length(args)>=1){
              eval(parse(text=paste("plot(0,0,bty='n',type='n',axes=FALSE,xlab='',ylab='',xlim=x.lim,ylim=y.lim,",paste(names(args),args,sep="=",collapse=","),")")))
            } else {
              plot(0,0,bty='n',type='n',axes=FALSE,xlab='',ylab='',xlim=x.lim,ylim=y.lim)
            }
            rect(c(1-freqs[1:show.barcode]),c((show.barcode-1):0),rep(1,show.barcode),c(show.barcode:1),col=ifelse(names(freqs)[1:show.barcode] %in% barcodes,"orange","black"))
            axis(2,labels=names(freqs)[1:show.barcode],at=seq(show.barcode-0.5,0.5,-1),las=2,tick=FALSE,line = -2,cex.lab=0.8)
            axis(3,labels=rev(seq(0,1,0.5)),at=seq(0,1,0.5))
            mtext("barcodes frequency",side=3,at=0.5,line=3)
            mtext(p.main,side=1,cex=par("cex.main"))
            
            ## return
            invisible(freqs)
          })

## for ShortReadQ
## think of the same for othe ShortRead obj and maybe a more generic function...
setMethod(
          f="barcodePlot",
          signature="ShortReadQ",
          definition=function(obj,barcodes=c(),type=c("independant","within"),barcode.length=6,show.barcode=20,...){
            return(
                   barcodePlot(sread(obj),barcodes=barcodes,type=type,barcode.length=barcode.length,show.barcode=show.barcode,...)
                   )
          })
