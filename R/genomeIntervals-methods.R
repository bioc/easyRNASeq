## convert a genome intervals into a RangedData
## TODO find a way to check and keep unexpected slots.
## or at least warn they would be ignored

setAs("Genome_intervals","RangedData",function(from){
  universe="intervals"
  
  # first check
  if (!is(from, "Genome_intervals")){stop("'from' must be a Genome_intervals object")}
  
  # create the ranges
  ranges<-IRanges(start=from[,1],end=from[,2])

  ## create the space
  ## drop the original levels
  space = as.character(seq_name(from))
  names(ranges) <- seq(along=ranges)
  
  # create the values
  slotToKeep<-na.omit(match(c("strand","exon","feature","intron","transcript","transcript.name","gene","gene.name"),names(from@annotation)))
  values <- as(data.frame(from@annotation[slotToKeep]),"DataFrame")
  elementMetadata(values)<-as(data.frame(labelDescription=names(from@annotation)[slotToKeep],row.names=names(from@annotation)[slotToKeep]),"DataFrame")

  # important before splitting
  rownames(values) <- names(ranges)
  
  # now split properly
  if (length(unique(space)) > 1) {
    ranges <- split(ranges, space)
    values <- split(values, space)
  } else {
    ranges <- RangesList(ranges)
    values <- SplitDataFrameList(values, compress = TRUE)
    names(ranges) <- unique(space)
    names(values) <- names(ranges)
  }
      
  if (!is.null(universe) && !isSingleString(universe))
    stop("'universe' must be a single string")
  universe(ranges) <- universe
  return(new(
      "RangedData",
      ranges = ranges,      
      values = values))
})

