## extend IRanges
setMethod(
          f="ranges",
          signature="RNAseq",
          definition=function(x){
            ranges(genomicAnnotation(x))
          })

