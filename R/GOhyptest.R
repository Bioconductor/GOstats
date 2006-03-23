##now we just sample K from the LLs and generate the GO graph

## GOHyperG moved to the Category package

GOKEGGHyperG <- function (x, lib = "hgu95av2", what = "MF", universe=NULL)
{
    .Defunct("geneKeggHyperGeoTest", "Category")

    getDataEnv <- function(name, lib) {
        get(paste(lib, name, sep = ""), mode = "environment")
    }



    require(lib, character.only = TRUE) || stop("need data package", lib)

    if (any(duplicated(x)))
      stop("input IDs must be unique")
    
    ################  MODIFIED #################
    match.arg(what, c("MF", "BP", "CC", "KEGG"))
    ############################################


    if( lib == "YEAST")
      cLLs <- unique(unlist(as.list(YEASTGO2ALLPROBES)))
    else
      cLLs <- unique(unlist(as.list(getDataEnv("LOCUSID", lib))))

    if(!is.null(universe)){
        if(any(duplicated(universe))){
            stop("universe must be unique")
        }
        goodUniverse <- universe %in% cLLs
        if(!all(goodUniverse)){
            ##FIXME:  any better statement for this?
            warning(paste("There are ", sum(!goodUniverse)), " out of ", length(goodUniverse), " elements in the universe which are not in the ", lib)
        }
        cLLs <- universe[goodUniverse]
    }

    goodX <- x %in% cLLs
    if(!all(goodX)){
        warning(paste("There are ", sum(!goodX), " out of ", length(goodX), " elememts in x which are not in the population"))
    }
    ourLLs <- x[x %in% cLLs]

    
    ################## MODIFIED ###########################
    if (what=="KEGG") {
        goV <- as.list(getDataEnv("PATH2PROBE", lib))
    } else {
        goV <- as.list(getDataEnv("GO2ALLPROBES", lib))
    }
    ######################################################

    ## use the subset of the whole GO stree
    if(!is.null(universe)){
        if( lib == "YEAST" ) {
            goV <- lapply(goV,
                          function(y){
                              if(all(unlist(y) %in% cLLs)) return(y)
                          })
            goV <- goV[which(sapply(goV, function(y) !is.null(y) ))]
        }else{
            goV <- lapply(goV,
                          function(y){
                              lls = unique(unlist(mget(y, getDataEnv("LOCUSID", lib),
                                ifnotfound=NA)))
                              if(all(lls %in% cLLs)) return(y)
                          })
            goV <- goV[which(sapply(goV, function(y) !is.null(y) ))]
        }
    }else{
        ## remove some non-interested samples in order to increase the speed
        if( lib == "YEAST" ){
            whWeHave <- sapply(goV,
                               function(y) {
                                   if( is.na(y) || length(y) == 0 )
                                     return(FALSE)
                                   any(ourLLs %in% y)
                               }
                               )
        }else{
            whWeHave <- sapply(goV,
                               function(y) {
                                   if( is.na(y) || length(y) == 0 )
                                     return(FALSE)
                                   lls = unique(unlist(mget(y, getDataEnv("LOCUSID", lib),
                                     ifnotfound=NA)))
                                   any(ourLLs %in% lls)
                               }
                               )
            
            goV <- goV[whWeHave]
        }
    }

    ##################  MODIFIED #########################
    if (what!="KEGG") {
        dropAll = match("all", names(goV), 0)
        if( dropAll > 0 )
          goV = goV[-dropAll]

        goCat <- unlist(getGOOntology(names(goV)))
        goodGO <- goCat == what
        mm <- match(names(goCat), names(goV))
        mm <- mm[goodGO]
        goV <- goV[mm]
    }
    ######################################################

    if( lib == "YEAST")
      goVcts = goV
    else
      goVcts = sapply(goV, function(x) {
          if(length(x) == 0 || is.na(x) ) return(NA)
          lls <- unique(unlist(mget(x, getDataEnv("LOCUSID", lib))))
          lls<-lls[!is.na(lls)]
          lls
      })

    bad <- sapply(goVcts, function(x) (length(x) == 1 && is.na(x)))
    goVcts = goVcts[!bad]
    goV = goV[!bad] ##we are going to output this...

    ##find out how many map which is m+n in phyper-speak
    cLLs <- unique(unlist(goVcts))
    nLL <- length(cLLs)

    goCounts <- sapply(goVcts, length)

    ours <- unique(ourLLs[!is.na(ourLLs)])
    ##need to get the number of interesting genes - these
    ##are the ones supplied that map to some GO term
    whGood <- ours[ours %in% cLLs]


    nInt = length(whGood)
    if( nInt == 0 )
       warning("no interesting genes found")

    useCts <- sapply(goVcts, function(x) sum(whGood %in% x))

    ##need the -1 because we are asking for evidence as extreme or
    ##more extreme than what we say - and we are using the upper tail
    pvs <- phyper(useCts-1, nInt, nLL-nInt, goCounts, lower.tail=FALSE)

    ord <- order(pvs)

        ################################  MODIFIED
##################################################
    if (what!="KEGG") {
        return(list(pvalues = pvs[ord], goCounts = goCounts[ord], chip = lib, go2Affy =
goV,
        intCounts = useCts[ord], numLL = nLL, numInt = nInt, intLLs = x))
    } else {
        return(list(pvalues = pvs[ord], keggCounts = goCounts[ord], chip = lib, kegg2Affy =
goV,
        intCounts = useCts[ord], numLL = nLL, numInt = nInt, intLLs = x))
    }
    ################################  MODIFIED
##################################################

}
