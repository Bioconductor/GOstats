##now we just sample K from the LLs and generate the GO graph

GOHyperG <- function(x, lib="hgu95av2", what="MF") {

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
    cLLs <- unlist(as.list(getDataEnv("LOCUSID", lib)))

    ##which ones do we have - don't worry about duplicates just yet
    ourLLs <- cLLs[match(x, cLLs)]

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

    whWeHave <- sapply(goV, function(y) {
        if( is.na(y) || length(y) == 0 )
            return(FALSE)
        lls = unique(unlist(mget(y, getDataEnv("LOCUSID", lib),
                        ifnotfound=NA)))
        any(x %in% lls) })

    goV <- goV[whWeHave]

    ##determine which of the GO terms are in the category we are
    ##working on
    goCat <- unlist(getGOCategory(names(goV)))
    goodGO <- goCat == what
    ##not every thing in goV has a category
    mm <- match(names(goCat), names(goV))
    mm <- mm[goodGO]
    goV <- goV[mm]


    ##this is the slow part
    goVcts = sapply(goV, function(x) {
        if(length(x) == 0 || is.na(x) ) return(NA)
        lls <- unique(unlist(mget(x, getDataEnv("LOCUSID", lib))))
        lls<-lls[!is.na(lls)]
        lls
    })
    bad <- sapply(goVcts, function(x) (length(x) == 1 && is.na(x)))
    goVcts = goVcts[!bad]

    ##find out how many map which is m+n in phyper-speak
    cLLs <- unique(unlist(goVcts))
    nLL <- length(cLLs)

    goCounts <- sapply(goVcts, length)


    ourLLs <- unique(ourLLs[!is.na(ourLLs)])
    ours <-ourLLs[!duplicated(ourLLs)]
    ##need to get the number of interesting genes - these
    ##are the ones supplied that map to some GO term
    whGood <- ours[ours %in% cLLs]


    nInt = length(whGood)
    if( nInt == 0 )
       warning("no interesting genes found")

    useCts <- sapply(goVcts, function(x) sum(whGood %in% x))

    ##need the -1 because we are asking for evidence as extreme or
    ##more extremen than what we say - and we are using the upper tail
    pvs <- phyper(useCts-1, nInt, nLL-nInt, goCounts, lower.tail=FALSE)

    ord <- order(pvs)
    return(list(pvalues=pvs[ord], goCounts=goCounts[ord],
                intCounts=useCts[ord], numLL=nLL, numInt = nInt))
}

    ##
