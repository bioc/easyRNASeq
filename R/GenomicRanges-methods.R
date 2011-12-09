## to extend GenomicRanges
setMethod(
          f="colnames",
          signature="GenomicRanges",
          definition=function(x, do.NULL = TRUE, prefix = "col"){
            if (length(x) == 0){
              return(character(0))
            } else {              
              return(switch(class(x),
                            "GRanges" = {
                              colnames(elementMetadata(x), do.NULL = do.NULL, prefix = prefix)
                            },
                            "GRangesList" = {
                              colnames(elementMetadata(x[[1]]), do.NULL = do.NULL, prefix = prefix)
                            }))
            }
          })

