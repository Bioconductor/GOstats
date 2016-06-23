##########################################
##
##  Copyright 2005 James W. MacDonald 
##
##  A function to output a table from running
##  GOHyperG() on a set of Affymetrix probe IDs
##
##
##  Modified 6-9-05 to output the table
##  Added hyperG2Affy 10-10-2005
##   - This function takes GOIDs form hyperGtable and outputs
##     a list giving the Affy IDs associated with each GOID
##
##
##  Moved hyperGtable and hyperG2Affy to GOstats 2-23-06
##
##  Added universe argument 3-22-06
##
## Made hyperGtable and hyperG3Affy defunct 10-2-07
###########################################

probeSetSummary <- function(result, pvalue, categorySize, sigProbesets,
                            ids = "ENTREZID") {
    if (!is(result, "GOHyperGResult"))
      stop("result must be a GOHyperGResult instance (or subclass)")
    ## build reverse map
    ps2egEnv <- getAnnMap(ids, annotation(result))
    psids <- ls(ps2egEnv)
    ## each psid maps to _exactly_ one egid
    ## Need to remove NAs.  This is where the database stuff
    ## will be _much_ easier.
    havePsids <- eapply(ps2egEnv, function(ids) {
        if (length(ids) == 1 && is.na(ids))
          FALSE
        else
          TRUE
    })
    ord <- order(names(havePsids))
    psids <- psids[unlist(havePsids)[ord]]
    eg2ps <- split(psids, unlist(mget(psids, ps2egEnv)))

    if (missing(pvalue))
      pvalue <- pvalueCutoff(result)
    if (missing(categorySize))
      categorySize <- NULL
    ## XXX FIXME:
    summary <- Category:::XXX_getSummaryGeneric_XXX()
    goids <- summary(result, pvalue, categorySize)[,1]
    sigegids <- geneIds(result)
    egids <- geneIdUniverse(result)[goids]
    psetids <- lapply(egids, function(ids) {
        ids <- as.character(ids)
        have <- ids %in% sigegids
        ans <- eg2ps[ids[have]]
        if (length(ans))
          unlist(ans)
        else
          NULL
    })
    psetidsNULL <- sapply(psetids, is.null)
    psetids <- psetids[!psetidsNULL]
    if(missing(sigProbesets)){
        selectedProbeSetIDs <- names(geneIds(result))
        if(is.null(selectedProbeSetIDs))
            warning(paste("The vector of geneIds used to create the GOHyperGParams",
                          "object was not a named vector.\nIf you want to know the",
                          "probesets that contributed to this result\neither use",
                          "a named vector for geneIds, or pass a vector of probeset IDs\n",
                          "via sigProbesets.\n", sep=""),
                    call.=FALSE)
    }else{
        selectedProbeSetIDs <- sigProbesets
    }
    selectedInd <- lapply(psetids, function(ids) {
        ids %in% selectedProbeSetIDs
    })
    egids <- lapply(psetids, function(ids) unlist(mget(ids, ps2egEnv)))
    out <- vector(mode="list", length=length(psetids))
    for (i in seq(along=psetids)) {
        out[[i]] <- data.frame(EntrezID=egids[[i]],
                               ProbeSetID=psetids[[i]],
                               selected=as.integer(selectedInd[[i]]),
                               row.names=NULL,
                               stringsAsFactors=FALSE)
    }
    names(out) <- names(psetids)
    out
}


setMethod("summary", signature(object="GOHyperGResult"),
          function(object, pvalue=pvalueCutoff(object),
                   categorySize=NULL, htmlLinks=FALSE) {
              AMIGO_URL <- "http://amigo.geneontology.org/amigo/term/%s"
              df <- callNextMethod(object=object, pvalue=pvalue,
                                   categorySize=categorySize)
              if (nrow(df) == 0)  {
                  df$Term <- character(0)
                  return(df)
              }
              goIds <- df[[1]]
              goTerms <- sapply(mget(goIds, GOenv("TERM")), Term)
              if (htmlLinks) {
                  goIdUrls <- sapply(goIds,
                                     function(x) sprintf(AMIGO_URL, x))
                  goTerms <- paste('<a href="', goIdUrls, '">', goTerms,
                                   '</a>', sep="")
              }
              df$Term <- goTerms
              df
          })



setMethod("htmlReport", signature(r="GOHyperGResult"),
          function(r, file="", append=FALSE, label="",
                   digits=3, summary.args=list(htmlLinks=TRUE))
          {
              ## We define a method here to pass htmlLinks=TRUE via
              ## ... as the default.

              ## NB: it "should" be enough to call callNextMethod() with
              ## no arguments, but it doesn't work.
              callNextMethod(r=r, file=file, append=append, label=label,
                             digits=digits, summary.args=summary.args)
          })


