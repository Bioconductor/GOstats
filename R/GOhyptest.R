##now we just sample K from the LLs and generate the GO graph

GOHyperG <- function(x, lib="hgu95av2", what="MF", universe=NULL) {

   ##a helper function
    getDataEnv <- function(name, lib) {
        get(paste(lib, name, sep=""), mode="environment")
    }

    ##first some checking

    require(lib, character.only=TRUE) || stop("need data package", lib)

    if( any(duplicated(x)) )
        stop("input IDs must be unique")

    match.arg(what, c("MF", "BP", "CC"))

    ##get the unique LocusLink IDs - GO is mapped to LL not to probe
    ##except for YEAST where it maps to common names

    if( lib == "YEAST")
      cLLs <- unique(unlist(as.list(YEASTGO2ALLPROBES)))
    else if (lib == "Universe")
      cLLs <- universe
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
    
    ##which ones do we have - don't worry about duplicates just yet
    goodX <- x %in% cLLs
    if(!all(goodX)){
        warning(paste("There are ", sum(!goodX), " out of ", length(goodX), " elememts in x which are not in the population"))
    }
    ourLLs <- x[x %in% cLLs]

    ##map to the unique LLs for each GO term
    goV <- as.list(getDataEnv("GO2ALLPROBES", lib))

    ##we need to reduce this, so first we only use the terms our
    ##symbols map to
    ##FIXME: this should work but the mappings are broken

    #a2GO <- mget(x, getDataEnv("GO", lib), ifnotfound=NA)
    #a2GOsz <- mget(unlist(sapply(a2GO, names)),
    #               env=hgu95av2GO2ALLPROBES, ifnotfound=NA)
    #aa <- sapply(a2GOsz, length)
    #a2GO <- unique(unlist(sapply(a2GO, names)))

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

    ##fix the all node problem
    dropAll = match("all", names(goV), 0)
    if( dropAll > 0 )
        goV = goV[-dropAll]

    ##Determine which of the GO terms are in the category we are
    ##working on
    ## fixed by Ting: getGOOntology will return a vector, not list
    ## goCat <- unlist(getGOOntology(names(goV)))
    goCat <- getGOOntology(names(goV))
    goodGO <- goCat == what
    ##not every thing in goV has a category
    mm <- match(names(goCat), names(goV))
    mm <- mm[goodGO]
    goV <- goV[mm]


    ##this is the slow part
    if( lib == "YEAST")
        goVcts = goV
    else
        goVcts = sapply(goV, function(y) {
            if(length(y) == 0 || is.na(y) ) return(NA)
            lls <- unique(unlist(mget(y, getDataEnv("LOCUSID", lib))))
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

    useCts <- sapply(goVcts, function(y) sum(whGood %in% y))

    ##need the -1 because we are asking for evidence as extreme or
    ##more extreme than what we say - and we are using the upper tail
    pvs <- phyper(useCts-1, nInt, nLL-nInt, goCounts, lower.tail=FALSE)

    ord <- order(pvs)
    return(list(pvalues=pvs[ord], goCounts=goCounts[ord], chip=lib,
                go2Affy=goV, intCounts=useCts[ord], numLL=nLL,
                numInt = nInt, intLLs=x))
}

    ##

GOKEGGHyperG <- function (x, lib = "hgu95av2", what = "MF", universe=NULL)
{

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
