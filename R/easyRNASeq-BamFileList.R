#' Get a BamFileList from a list of filenames
#'
#' A utility function to create a \code{linkS4class{BamFileList-class}}{BamFileList}
#' object from a set of filenames. The filenames need to contain the file path if they
#' are not in the working directory.

#' @name getBamFileList
#' @rdname easyRNASeq-BamFileList
#' @aliases getBamFileList,character,character-method
#' getBamFileList,character,missing-method
#' @param filenames a character vector containing fully defined BAM file filenames
#' @param indexnames a character vector containing fully defined BAM index file filenames
#' @return a \code{linkS4class{BamFileList-class}}{BamFileList}
#' @seealso \code{linkS4class{BamFileList-class}}{BamFileList}
#' \code{\link[base:list.files]{dir}}
#' @examples
#' # tutorial data - store the data in the BiocCache
#' tdir <- tutorialData()
#'
#'   # creating a BamFileList using a directory and pattern
#'   # using filenames (from the Bioc cache)
#'   filenames <- dir(tdir,pattern="[A,C,T,G]{6}\\.bam$",full.names=TRUE)
#'   indexnames <- sapply(paste0(sub(".*_","",basename(filenames)),".bai"),fetchData)
#'   bfl <- getBamFileList(filenames,indexnames)
#'
#'   # get them recursively
#'   filenames <- dir(path=tdir,pattern="[A,C,T,G]{6}\\.bam$",
#'                    full.names=TRUE,recursive=TRUE)
#'   indexnames <- sapply(paste0(sub(".*_","",basename(filenames)),".bai"),fetchData)
#'   bfl <- getBamFileList(filenames,indexnames)
#'
setMethod(f="getBamFileList",
          signature=c("character","character"),
          definition=function(filenames=character(0),
                              indexnames=character(0)){

            # check the files
            if(length(filenames) == 0){
              stop("You need to provide a character vector of 'filenames'. ")
            }
              if(length(indexnames) == 0){
                  stop("You need to provide a character vector of 'indexnames'. ")
              }

            # check that the file(s) exist(s)
            sel <- file.exists(filenames)
            if(any(!sel)){
              stop(paste("The file(s)",
                         paste(filenames[!sel],collapse = " and ")),
                   " do(es) not exist.")
            }

            sel <- file.exists(indexnames)
            if(any(!sel)){
                stop(paste("The index file(s)",
                           paste(indexnames[!sel],collapse = " and ")),
                     " do(es) not exist.")
            }

            # check that the bai are not part of the filenames
            if(length(grep("\\.bam\\.bai$",filenames))>0){
                warning(paste("You either have provided BAM index files (.bai) as",
                              "part of the filenames argument or your pattern matched",
                              "BAM index files! Removing these from the files",
                              "list to process."))
                filenames <- filenames[!grepl("\\.bam\\.bai$",filenames)]
            }

            # check if we are having any bam files
            if(length(filenames)==0){
              stop("There are no file to process.")
            }

            # create the BamFileList
            bfl <- BamFileList(filenames,index=indexnames)
            names(bfl) <- basename(filenames)

            # validate them
            validate(bfl)

            # return
            return(bfl)
          })

setMethod(f="getBamFileList",
          signature=c("character","missing"),
          definition=function(filenames=character(0)){

              # check if we have index with bai
              indexes <- paste(filenames,"bai",sep=".")
              sel <- file.exists(indexes)
              if(any(!sel)){
                  stop(paste("Index files (bai) are required. They are missing for the files: ",
                             paste(filenames[!sel],collapse = " and "),
                             ". Use the Rsamtools indexBam function or the samtools index command line utility to create them.",sep=""))
              }
              getBamFileList(filenames,indexes)
          })

