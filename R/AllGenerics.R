setGeneric("getGoGraph", function(p, goIds) standardGeneric("getGoGraph"))

setGeneric("goDag", function(r) standardGeneric("goDag"))

setGeneric("geneIdUniverse", function(r) standardGeneric("geneIdUniverse"))

setGeneric("condGeneIdUniverse",
           function(r) standardGeneric("condGeneIdUniverse"))


## Accessors for GOHyperGResult

setGeneric("htmlReport", function(r, file="", append=TRUE, label="", ...)
           standardGeneric("htmlReport"),
           signature=c("r"))

## Create generics for non-generics defined in base
setGeneric("summary")
