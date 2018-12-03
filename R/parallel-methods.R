##' parallel additional methods
##'
##' Functions defined in the easyRNASeq package that enhance the parallel package.
##'
##' The parallelize function ease the use of the parallel package. If the number
##' of nodes provided by the user is 1, then a simple 'lapply' is used, otherwise
##' a cluster object is created and the object dispatched for parallelization.
##'
##' @rdname parallel-methods
##' @name parallel additional methods
##' @aliases parallelize parallelize,list,function-method
##' parallelize,vector,function-method
##' parallelize,GRangesList,function-method
##' parallelize,BamFileList,function-method
##' @param fun the function to be applied in parallel
##' @param nnodes the number of nodes to use
##' @param obj the object which processing has to be parallelizes
##' @param ... additional arguments passed to the function \code{fun}
##' @return the result of the \code{\link[parallel:clusterApply]{clusterApply}} function.
##' @author Nicolas Delhomme
##' @seealso \code{\link[parallel:clusterApply]{clusterApply}}
##' \code{makePSOCKcluster} and \code{stopCluster} in \code{\link[parallel]{makeCluster}}
##' @keywords methods
##' @examples
##' parallelize(list(a<-c(1,2),b<-c(2,1)),sum,nnodes=1)
setMethod(f = "parallelize",
          signature = c("list","function"),
          definition = function(obj,fun,nnodes=1,...){
            .parallelize(obj,fun,nnodes,...)
          })

setMethod(f = "parallelize",
          signature = c("vector","function"),
          definition = function(obj,fun,nnodes=1,...){
            obj <- as.list(obj)
            names(obj) <- unlist(obj)
            parallelize(obj,fun,nnodes,...)
          })

setMethod(f = "parallelize",
          signature = c("GRangesList","function"),
          definition = function(obj,fun,nnodes=1,...){
            .parallelize(obj,fun,nnodes,...)
          })

setMethod(f = "parallelize",
          signature = c("BamFileList","function"),
          definition = function(obj,fun,nnodes=1,...){
            .parallelize(obj,fun,nnodes,...)
          })

".parallelize" <- function(obj,fun,nnodes=1,...){
  if(nnodes == 1){
    res <- lapply(obj,fun,...)
  } else {
    ## parallel is part of R, we just ensure it's loaded
    ## require("parallel")
    ## save the names
    nams <- names(obj)
    cluster <- makePSOCKcluster(nnodes)
    res <- clusterApply(cluster,obj,fun,...)
    stopCluster(cluster)
    names(res) <- nams
  }
  return(res)
}

