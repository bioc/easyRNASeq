## show, just call print.
setMethod(
          f="show",
          signature="RNAseq",          
          definition=function(object){
            print(object)
          })
