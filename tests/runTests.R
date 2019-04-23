# get the example data
library(easyRNASeq)
tutorialData()

# set the env.var
#TUTORIAL.DATA <- get("TUTORIAL.DATA",envir=as.environment("package:easyRNASeq"))

# run the tests
BiocGenerics:::testPackage("easyRNASeq")

# cleanup
# removebfc(easyRNASeq:::.get_cache(),ask=FALSE)
