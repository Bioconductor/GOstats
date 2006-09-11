setMethod("goHyperGeoTest",
          signature(p="GeneGoHyperGeoTestParams"), 
          function(p) {
              p@universeGeneIds <- universeBuilder(p)
              selected <- intersect(p@geneIds, p@universeGeneIds)
              p@geneIds <- selected
              cat2Entrez <- categoryToEntrezBuilder(p)
              ## build the GO graph for the relevant GO nodes
              goIds <- names(cat2Entrez)
              goDag <- getGoGraph(p, goIds)
              nodeDataDefaults(goDag, "pvalue") <- 1
              nodeDataDefaults(goDag, "geneIds") <- numeric(0)
              nodeDataDefaults(goDag, "condGeneIds") <- numeric(0)
              nodeDataDefaults(goDag, "oddsRatio") <- numeric(0)
              nodeDataDefaults(goDag, "expCount") <- numeric(0)
              ## store the Entrez Gene IDs as attrs on the GO DAG.
              nodeData(goDag, n=names(cat2Entrez),
                       attr="geneIds") <- cat2Entrez
              ## now iterate leaves first doing tests and conditioning
              ## on all significant children.
              ## FIXME: consider replacing with RBGL tsort?
              needsProc <- goDag
              complete <- character(0)
              SIGNIF <- p@pvalueCutoff
              while (length(nodes(needsProc))) {
                  numKids <- sapply(edges(needsProc), length)
                  noKids <- names(numKids[numKids == 0])
                  curCat2Entrez <- cat2Entrez[noKids]
                  if (p@conditional) {
                      curCatKids <- edges(goDag)[names(curCat2Entrez)]
                      curCatKids <- removeLengthZero(curCatKids)
                      if (length(curCatKids)) { ## sanity check
                          ## they should be all complete
                          stopifnot(all(unlist(curCatKids) %in% complete))
                      }
                      curCat2Entrez <- removeSigKidGenes(curCatKids, goDag,
                                                         curCat2Entrez,
                                                         SIGNIF, cat2Entrez)
                      ## Store the conditioned cat => entrez map
                      nodeData(goDag, n=names(curCat2Entrez),
                               attr="condGeneIds") <- curCat2Entrez
                  }
                  stats <- getHyperGeoPvalues(p, curCat2Entrez,
                                              cat2Entrez, selected)
                  ##odds_ratio <- getOddsRatio(p, curCat2Entrez, selected)

                  ## store the pvals, mark these nodes as complete,
                  ## then compute the next set of nodes to do.  
                  noKids <- names(curCat2Entrez)
                  ## drop names on pvals to avoid weird names upon unlisting
                  nodeData(goDag, n=noKids,
                           attr="pvalue") <- as.numeric(stats$p)
                  nodeData(goDag, n=noKids,
                           attr="oddsRatio") <- as.numeric(stats$odds)
                  nodeData(goDag, n=noKids,
                           attr="expCount") <- as.numeric(stats$expected)
                  complete <- c(complete, noKids)
                  hasKids <- names(numKids[numKids > 0])
                  needsProc <- subGraph(hasKids, needsProc)
              } ## end while
              pvals <- unlist(nodeData(goDag, attr="pvalue"))
              new("GeneGoHyperGeoTestResult",
                  goDag=goDag,
                  annotation=p@annotation,
                  geneIds=p@geneIds,
                  testName=categoryName(p),
                  testDirection=p@testDirection,
                  pvalueCutoff=p@pvalueCutoff,
                  pvalue.order=order(pvals))
          })


removeLengthZero <- function(x) {
    wanted <- sapply(x, function(z) length(z) > 0)
    x[wanted]
}


setMethod("getGoGraph", signature(p="GeneGoHyperGeoTestParams",
                                  goIds="character"),
          function(p, goIds) {
              ## FIXME: ':::'
              goEnv <- Category:::getDataEnv(paste(p@ontology,
                                                   "CHILDREN", sep=""),
                                             "GO")
              gobpkids <- eapply(goEnv, function(x) {
                  intersect(x, goIds)
              })
              gobpkids <- l2e(gobpkids)
              GOGraph(goIds, gobpkids)
          })


removeSigKidGenes <- function(curCatKids, goDag, curCat2Entrez, SIGNIF,
                              cat2Entrez) {
    if (length(curCatKids)) {
        ## keep only those kids with SIGNIF pvalue
        curCatKids <- lapply(curCatKids, function(x) {
            pvKids <- nodeData(goDag, n=x, attr="pvalue")
            idx <- which(pvKids < SIGNIF)
            if (length(idx))
              x[idx]
            else
              character(0)
        })
        curCat2EntrezCond <- list()
        for (goid in names(curCat2Entrez)) {
            ## remove entrez ids that came from 
            ## SIGNIF children
            kids <- curCatKids[[goid]]
            if (length(kids)) {
                kidEgs <- unlist(cat2Entrez[kids])
                newEgs <- setdiff(curCat2Entrez[[goid]], kidEgs)
                ## newEgs may be length 0
                curCat2EntrezCond[[goid]] <- newEgs
            } else {
                curCat2EntrezCond[[goid]] <- curCat2Entrez[[goid]]
            }
        }
        curCat2Entrez <- curCat2EntrezCond
    }
    curCat2Entrez
}


getHyperGeoPvalues <- function(p, curCat2Entrez, cat2Entrez, selected) {
    ## Here is how we conceptualize the test:
    ##
    ## The urn contains genes from the gene universe.  For each GO ID,
    ## genes annotated at the given GO term are white and the rest
    ## black.
    ##
    ## The number drawn is the size of the selected gene list.  The
    ## number of white drawn is the size of the intersection of the
    ## selected list and the GO list.
    ##
    ## In the conditional case, the GO ID annotation set has been
    ## reduced and we also adjust the selected list (num drawn) and gene
    ## universe according to what was removed by the conditioning.
    ##
    ##          inGO    notGO
    ##          White   Black
    ## selected  n11     n12
    ## not       n21     n22
    ##
    
    if (p@conditional) {
        cat2RemovedEntrez <- lapply(names(curCat2Entrez),
                                    function(goid) {
                                        setdiff(cat2Entrez[[goid]],
                                                curCat2Entrez[[goid]])
                                    })

        ## White balls removed from urn by conditioning
        numSelectedRemoved <- sapply(cat2RemovedEntrez,
                                     function(x) sum(selected %in% x))

        ## Black balls removed from urn by conditioning
        numOtherRemoved <- listLen(cat2RemovedEntrez) - numSelectedRemoved
    } else {
        numSelectedRemoved <- rep(0, length(curCat2Entrez))
        numOtherRemoved <- numSelectedRemoved
    }

    ## Num white drawn (n11)
    numWdrawn <- sapply(curCat2Entrez, 
                        function(x) sum(selected %in% x))

    ## Num white
    numW <- listLen(curCat2Entrez)
    
    ## Num black
    numB <- length(p@universeGeneIds) - numOtherRemoved - numW

    ## Num drawn
    numDrawn <- length(selected) - numSelectedRemoved

    n21 <- numW - numWdrawn
    n12 <- numDrawn - numWdrawn
    n22 <- numB - n12

    odds_ratio <-  (numWdrawn * n22) / (n12 * n21)

    expected <- (numWdrawn + n12) * (numWdrawn + n21)
    expected <- expected / (numWdrawn + n21 + n21 + n22)

    if (p@testDirection == "over") {
        ## take the -1 because we want evidence for as extreme or more
        pvals <- phyper(numWdrawn - 1, numW, numB,
                        numDrawn, lower.tail=FALSE)
    } else {
        pvals <- phyper(numWdrawn, numW, numB,
                        numDrawn, lower.tail=TRUE)
    }
    list(p=pvals, odds=odds_ratio, expected=expected)
}


GOHG <- function(entrezGeneIds, lib, ontology, universe=NULL,
                 conditional=FALSE, testDirection="over",
                 pvalueCutoff=0.01)
{
    if (missing(universe) || is.null(universe))
      universe <- character(0)

    if (conditional && missing(pvalueCutoff))
      stop("conditional computation requires ",
           sQuote("pvalueCutoff"), " to be specified")
    if (missing(pvalueCutoff))
      pvalueCutoff <-  0.01
    
    params <- new("GeneGoHyperGeoTestParams",
                  geneIds=entrezGeneIds,
                  universeGeneIds=universe,
                  annotation=lib,
                  ontology=ontology,
                  testDirection=testDirection,
                  pvalueCutoff=pvalueCutoff)
    
    goHyperGeoTest(params)
}
