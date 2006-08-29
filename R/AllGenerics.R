setGeneric("getGoGraph", function(p, goIds) standardGeneric("getGoGraph"))

setGeneric("goDag", function(r) standardGeneric("goDag"))

setGeneric("goHyperGeoTest", 
           function(p) {
               standardGeneric("goHyperGeoTest")
           }, valueClass="GeneGoHyperGeoTestResult")
