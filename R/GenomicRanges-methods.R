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
                              colnames(x, do.NULL = do.NULL, prefix = prefix)
                            },
                            {stop(paste("Unknown class:",class(x)))}))
            }
          })

## R 2.15.0 introduced the GRangesList signature
## TODO not sure if the switch in the method above is still required
setMethod(
          f="colnames",
          signature="GRangesList",
          definition=function(x, do.NULL = TRUE, prefix = "col"){
            if (length(x) == 0){
              return(character(0))
            } else {              
              colnames(elementMetadata(x[[1]]), do.NULL = do.NULL, prefix = prefix)              
            }
          })

