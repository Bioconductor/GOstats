##Copyright R. Gentleman, 2004
##simple functions to get Evidence codes

##get then GO term names for a particular (sub)ontology
getOntology = function(inlist, ontology=c("MF", "BP", "CC")) {
   which = match.arg(ontology)
   lapply(inlist, function(x) {
       onts = sapply(x, function(z) z$Ontology)
       onts = onts[!is.na(onts)]
       unique(names(x[onts %in% which]))
   })}


##get GO evidence codes
getEvidence = function(inlist) {
    lapply(inlist, function(x)
          sapply(x, function(z) z$Evidence))
}

##drop a specified set of evidence codes
dropECode = function(inlist, code = "IEA") {
    lapply(inlist, function(x) {
        hasCode = sapply(x, function(z) z$Evidence)
        badVals = hasCode %in% code
        x[!badVals]
    })
 }
