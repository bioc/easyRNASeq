## test some print function(s)
### =========================
## AnnotParam related
### =========================
"test_printAnnotParamObject" <- function(){

  obj <- new("AnnotParamObject")

  ## the printed value
  checkTrue(capture.output(print(obj))=="Annotation provided manually ")
}
