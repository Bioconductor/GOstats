setGeneric("getGoGraph", function(p, goIds) standardGeneric("getGoGraph"))

setGeneric("goDag", function(r) standardGeneric("goDag"))

setGeneric("condHyperGeoTest", 
           function(p) {
               standardGeneric("condHyperGeoTest")
           }, valueClass="GeneGoHyperGeoTestResult")
