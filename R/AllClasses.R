## Not sure we need a prototype later on

## TODO should we have one readAnnotation slot
## containing all the genomicAnnoation, geneModel, and readIslands info?
## and have it more dynamic?

setClass(
         Class="RNAseq",
         representation=representation(
           chrSize="list",
           fileName="character",
           geneModel="RangedData",
           genomicAnnotation="Vector",
           librarySize="numeric",
           organismName="character",
           readCounts="list",
           readCoverage="RleList",
           readIslands="RangedData",
           readLength="integer"
           ),
         prototype=prototype(
           chrSize=list(),
           fileName=character(0),
           geneModel=RangedData(),           
           genomicAnnotation=RangedData(),
           librarySize=numeric(0),
           organismName=character(1),
           readCounts=list(),
           readCoverage=RleList(),           
           readIslands=RangedData(),
           readLength=36L
           )
##         ,
##          validity=function(obj){
##          }
         )

