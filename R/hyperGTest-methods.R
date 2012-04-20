setMethod("hyperGTest", signature(p="GOHyperGParams"),
          function(p) .hyperGTestInternal(p) )

.hyperGTestInternal <- function(p) {
  p <- makeValidParams(p)
  p@universeGeneIds <- universeBuilder(p)
  ## preserve names on geneIds
  p@geneIds <- p@geneIds[p@geneIds %in% p@universeGeneIds]
  ##FIXME: since the unique values in this list are the
  ## the universe, it might be nice to have that somewhere
  cat2Entrez <- categoryToEntrezBuilder(p)
  ## build the GO graph for the relevant GO nodes
  goIds <- names(cat2Entrez)
  goDag <- getGoToChildGraph(p, goIds)
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
    stats <- .doHyperGTest(p, curCat2Entrez, cat2Entrez,
                           p@geneIds)
    
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
  ##modded for future generic GenSet baseed objects...
  if(class(p@datPkg)!="GeneSetCollectionDatPkg"){
    annotation = annotation(p)
  }else if(class(p@datPkg)=="GeneSetCollectionDatPkg"){
    annotation = "Based on a GeneSetCollection Object"
  }
  new("GOHyperGResult",
      goDag=goDag,
      annotation=annotation,
      geneIds=p@geneIds,
      testName=categoryName(p),
      testDirection=p@testDirection,
      pvalueCutoff=p@pvalueCutoff,
      pvalue.order=order(pvals),
      conditional=p@conditional)
}





removeLengthZero <- function(x) {
    wanted <- sapply(x, function(z) length(z) > 0)
    x[wanted]
}


getGoToChildGraph <- function(p, goIds) {
    ## We might want a method here so that we can dispatch on
    ## a SQLite-based GO package.
    goEnv <- GOenv(paste(p@ontology, "CHILDREN", sep=""))
    gobpkids <- eapply(goEnv, function(x) {
        intersect(x, goIds)
    })
    gobpkids <- list2env(gobpkids)
    GOGraph(goIds, gobpkids)
}


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
    
    params <- new("GOHyperGParams",
                  geneIds=entrezGeneIds,
                  universeGeneIds=universe,
                  annotation=lib,
                  ontology=ontology,
                  testDirection=testDirection,
                  pvalueCutoff=pvalueCutoff)
    
    hyperGTest(params)
}
