##Copyright R. Gentleman, 2004
##simple functions to get Evidence codes

##get then GO term names for a particular (sub)ontology
getOntology = function(inlist, ontology=c("MF", "BP", "CC")) {
   which = match.arg(ontology)
   onts = sapply(inlist, function(z) z$Ontology)
   onts = onts[!is.na(onts)]
   unique(names(inlist[onts %in% which]))
}


##get GO evidence codes
getEvidence = function(inlist) 
     sapply(inlist, function(z) z$Evidence)

##drop a specified set of evidence codes
dropECode = function(inlist, code = "IEA") {
    hasCode = sapply(inlist, function(z) z$Evidence)
    badVals = hasCode %in% code
    inlist[!badVals]
}
