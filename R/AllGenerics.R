setGeneric("getGoGraph", function(p, goIds) standardGeneric("getGoGraph"))

setGeneric("goDag", function(r) standardGeneric("goDag"))

setGeneric("geneIdUniverse", function(r) standardGeneric("geneIdUniverse"))

setGeneric("condGeneIdUniverse",
           function(r) standardGeneric("condGeneIdUniverse"))

setGeneric("goHyperGeoTest", 
           function(p) {
               standardGeneric("goHyperGeoTest")
           }, valueClass="GOHyperGResult")

setGeneric("oddsRatios", function(r) standardGeneric("oddsRatios"))

setGeneric("expectedCounts",
           function(r) standardGeneric("expectedCounts"))

## Accessors for GOHyperGResult

setGeneric("isConditional", function(r) standardGeneric("isConditional"))

setGeneric("htmlReport", function(r, file="", append=TRUE, label="", ...)
           standardGeneric("htmlReport"),
           signature=c("r"))

## Create generics for non-generics defined in base
setGeneric("summary")
