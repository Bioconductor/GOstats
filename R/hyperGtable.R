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
###########################################

hyperGtable <- function(probids, lib, type="MF", pvalue=0.05,
                        min.count=10, save = FALSE, output = TRUE,
                        filename = NULL, universe = NULL){
    .Deprecated("summary")
    message("use summary(result) where result is the output of hyperGTest")
  require(lib, quietly = TRUE, character.only = TRUE) || stop(paste("The ", lib, " package is required"))
  lls <- getLL(probids, lib)
  lls <- unique(lls)
  lls <- lls[!is.na(lls)]
  tmp <- GOHyperG(lls, lib, type, universe)
  index <- tmp$pvalues < pvalue & tmp$goCounts > min.count
  if(sum(index) == 0)
    stop(paste("There are no significant GO terms using a p-value of ",pvalue,
               " and a minimum count of ",  min.count,".\n", sep = ""), call. = FALSE)
  wh <- mget(names(tmp$pvalues[index]), GOTERM)
  tmp.terms <- sapply(wh, Term)
  out <- data.frame(names(tmp$pvalues[index]), tmp.terms, round(tmp$pvalues[index], 3), 
                       tmp$intCounts[index], tmp$goCounts[index])
  names(out) <- c("GO Term","Description","p-values",
                     "Number significant","Number on chip")
  if(output){
    if(!is.null(filename)){
      write.table(out, paste(filename, ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
    }else{
      stop("You need to give a filename in order to output this table!\n")
    }
  }
  if(save)
    out
}

hyperG2Affy <- function(probids, lib, type="MF", pvalue=0.05,
                        min.count=10, universe = NULL){
    .Deprecated("probeSetSummary")
  require(lib, quietly = TRUE, character.only = TRUE) || stop(paste("The ", lib, " package is required"))
  lls <- getLL(probids, lib)
  lls <- unique(lls)
  lls <- lls[!is.na(lls)]
  tmp <- GOHyperG(lls, lib, type, universe)
  index <- tmp$pvalues < pvalue & tmp$goCounts > min.count
  if(sum(index) == 0)
    stop(paste("There are no significant GO terms using a p-value of ",pvalue,
               " and a minimum count of ",  min.count,".\n", sep = ""), call. = FALSE)
  wh <- mget(names(tmp$pvalues[index]), GOTERM)
  tmp.terms <- sapply(wh, Term)
  index2 <- match(names(tmp.terms), names(tmp$go2Affy))
  out <- lapply(tmp$go2Affy[index2], function(x) x[x %in% probids])
  names(out) <- names(tmp.terms)
  out
}


probeSetSummary <- function(result, pvalue, categorySize) {
    if (!is(result, "GOHyperGResult"))
      stop("result must be a GOHyperGResult instance (or subclass)")
    ## build reverse map
    ps2egEnv <- Category:::getDataEnv("ENTREZID", annotation(result))
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
    goids <- summary(result, pvalue, categorySize)[,1]
    ## XXX: these are unconditional, not sure if we want the
    ##      condGeneIdUniverse here if the calculation used
    ##      the conditional calculation.
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
    selectedProbeSetIDs <- names(geneIds(result))
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
                   categorySize=NULL, htmlLinks=TRUE) {
              AMIGO_URL <- "http://www.godatabase.org/cgi-bin/amigo/go.cgi?view=details&search_constraint=terms&depth=0&query=%s"
              df <- callNextMethod(object=object, pvalue=pvalue,
                                   categorySize=categorySize)
              if (nrow(df) == 0)  {
                  df$Term <- character(0)
                  return(df)
              }
              goIds <- df[[1]]
              goTerms <- sapply(mget(goIds, GOTERM), Term)
              if (htmlLinks) {
                  goIdUrls <- sapply(goIds,
                                     function(x) sprintf(AMIGO_URL, x))
                  goTerms <- paste('<a href="', goIdUrls, '">', goTerms,
                                   '</a>', sep="")
              }
              df$Term <- goTerms
              df
          })
