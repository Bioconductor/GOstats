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
              ## store the Entrez Gene IDs as attrs on the GO DAG.
              nodeData(goDag, n=names(cat2Entrez),
                       attr="geneIds") <- cat2Entrez
              ## now iterate leaves first doing tests and conditioning
              ## on all significant children.
              ## FIXME: consider replacing with RBGL tsort?
              needsProc <- goDag
              complete <- character(0)
              SIGNIF <- p@pvalue.cutoff
              while (length(nodes(needsProc))) {
                  numKids <- sapply(edges(needsProc), length)
                  noKids <- names(numKids[numKids == 0])
                  curCat2Entrez <- cat2Entrez[noKids]
                  if (p@conditional) {
                      curCatKids <- edges(goDag)[names(curCat2Entrez)]
                      curCatKids <- removeLengthZero(curCatKids)
                  } else {
                      curCatKids <- character(0)
                  }
                  tmp <- removeSigKidGenes(curCatKids, goDag,
                                           complete, curCat2Entrez,
                                           SIGNIF, cat2Entrez)
                  curCat2Entrez <- tmp$curCat2Entrez
                  complete <- tmp$complete
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
                  test.direction=p@test.direction,
                  pvalue.cutoff=p@pvalue.cutoff,
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


removeSigKidGenes <- function(curCatKids, goDag, complete,
                              curCat2Entrez, SIGNIF, cat2Entrez) {
    if (length(curCatKids)) {
        ## they should be all complete (sanity check)
        stopifnot(all(unlist(curCatKids) %in% complete))
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
                if (length(newEgs))
                  curCat2EntrezCond[[goid]] <- newEgs
                else ## we skip the test, but mark complete
                  complete <- c(complete, goid)
            } else {
                curCat2EntrezCond[[goid]] <- curCat2Entrez[[goid]]
            }
        }
        curCat2Entrez <- curCat2EntrezCond
    }
    list(curCat2Entrez=curCat2Entrez, complete=complete)
}


getHyperGeoPvalues <- function(p, curCat2Entrez, selected) {
    numFound <- sapply(curCat2Entrez, 
                       function(x) sum(selected %in% x))
    numDrawn <- length(selected)
    ## num white in urn
    numAtCat <- sapply(curCat2Entrez, length)
    ## num black in urn
    numNotAtCat <- length(p@universeGeneIds) - numAtCat
    ## take the -1 because we want evidence for as extreme or
    ## more.
    ## NOTE: change to pyhper(numFound, ...) to
    ## compute underrepresentation.
    overUnder <- switch(p@test.direction,
                        over=1,
                        under=0,
                        stop("Unknown test.direction",
                             p@test.direction))
    pvals <- phyper(numFound-overUnder, numAtCat, numNotAtCat,
                    numDrawn, lower.tail=FALSE)
    pvals
}


GOHG <- function(entrezGeneIds, lib, ontology, universe=NULL,
                 conditional=FALSE, test.direction="over",
                 pvalue.cutoff=0.01)
{
    if (missing(universe) || is.null(universe))
      universe <- character(0)

    if (conditional && missing(pvalue.cutoff))
      stop("conditional computation requires ",
           sQuote("pvalue.cutoff"), " to be specified")
    if (missing(pvalue.cutoff))
      pvalue.cutoff <-  0.01
    
    params <- new("GeneGoHyperGeoTestParams",
                  geneIds=entrezGeneIds,
                  universeGeneIds=universe,
                  annotation=lib,
                  ontology=ontology,
                  test.direction=test.direction,
                  pvalue.cutoff=pvalue.cutoff)
    
    goHyperGeoTest(params)
}
