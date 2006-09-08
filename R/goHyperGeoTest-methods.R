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
                  pvals <- getHyperGeoPvalues(p, curCat2Entrez, selected)

                  ## store the pvals, mark these nodes as complete,
                  ## then compute the next set of nodes to do.  
                  noKids <- names(curCat2Entrez)
                  ## drop names on pvals to avoid weird names upon unlisting
                  nodeData(goDag, n=noKids, attr="pvalue") <- as.numeric(pvals)
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


getHyperGeoPvalues <- function(p, curCat2Entrez, selected) {
    numFound <- sapply(curCat2Entrez, 
                       function(x) sum(selected %in% x))
    numDrawn <- length(selected)
    ## num white in urn
    numAtCat <- sapply(curCat2Entrez, length)
    ## num black in urn
    numNotAtCat <- length(p@universeGeneIds) - numAtCat
    if (p@testDirection == "over") {
        ## take the -1 because we want evidence for as extreme or more
        pvals <- phyper(numFound - 1, numAtCat, numNotAtCat,
                        numDrawn, lower.tail=FALSE)
    } else {
        pvals <- phyper(numFound, numAtCat, numNotAtCat,
                        numDrawn, lower.tail=TRUE)
    }
    pvals
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
