setMethod("hyperGTest", "GOHyperGParams", function(p) {
    .hyperGTestInternal(p, "GOHyperGResult", getGoToChildGraph)
})

setMethod("hyperGTest", "OBOHyperGParams", function(p) {
    .hyperGTestInternal(p, "OBOHyperGResult", getGoToChildOBOgraph)
})

.hyperGTestInternal <- function(p, resultsClassName, childGraphFun) {
  p <- makeValidParams(p)
  p@universeGeneIds <- universeBuilder(p)
  ## preserve names on geneIds
  p@geneIds <- p@geneIds[p@geneIds %in% p@universeGeneIds]
  ##FIXME: since the unique values in this list are the
  ## the universe, it might be nice to have that somewhere
  cat2Entrez <- categoryToEntrezBuilder(p)
  ## build the GO graph for the relevant GO nodes
  goIds <- names(cat2Entrez)
  goDag <- childGraphFun(p, goIds)
  nodeDataDefaults(goDag, "pvalue") <- 1
  nodeDataDefaults(goDag, "geneIds") <- numeric(0)
  nodeDataDefaults(goDag, "condGeneIds") <- numeric(0)
  nodeDataDefaults(goDag, "oddsRatio") <- 1
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
  SIGODD <- p@orCutoff
  MINSIGSZE <- p@minSizeCutoff
  MAXSIGSZE <- p@maxSizeCutoff
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
                                         SIGNIF, SIGODD, MINSIGSZE, MAXSIGSZE,
                                         cat2Entrez)
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
  objparams <- list(Class=resultsClassName,
                    goDag=goDag,
                    annotation=annotation,
                    geneIds=p@geneIds,
                    testName=categoryName(p),
                    testDirection=p@testDirection,
                    pvalueCutoff=p@pvalueCutoff,
                    pvalue.order=order(pvals),
                    conditional=p@conditional)
  if (resultsClassName == "OBOHyperGResult") {
      desc <- sapply(p@datPkg@geneSetCollection[nodes(goDag)], description)
      objparams$gscDescriptions <- setNames(desc, nodes(goDag))
  }
  do.call("new", objparams)
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

## this is for OBO
getGoToChildOBOgraph <- function(p, goIds) {
  childrenEdges <- edges(as(t(as(p@datPkg@oboGraph, "matrix")), "graphNEL"))
  childrenEdges <- lapply(childrenEdges, intersect, goIds)
  childrenEdges <- list2env(childrenEdges)
  GOGraph(goIds, childrenEdges)
}


## ROBERT 30.7.2017. Added arguments SIGODD, MINSIGSZE, MAXSIGSZE
## to consider odds ratio and gene set size (minimum and maximum)
## when removing children from significant terms. Thus the notion
## of a "significant" term is extended to those meeting also these
## parameters, when given.
removeSigKidGenes <- function(curCatKids, goDag, curCat2Entrez, SIGNIF,
                              SIGODD, MINSIGSZE, MAXSIGSZE, cat2Entrez) {
    if (length(curCatKids)) {
        ## keep only those kids with SIGNIF pvalue
        ## SIGODD odds ratio, MINSIGSZE minimum gene set size
        ## and MAXSIGSZE maximum gene set size
        curCatKids <- lapply(curCatKids, function(x) {
            pvKids <- nodeData(goDag, n=x, attr="pvalue")
            orKids <- nodeData(goDag, n=x, attr="oddsRatio")
            szeKids <- lengths(nodeData(goDag, n=x, attr="geneIds"))
            idx <- which(pvKids < SIGNIF & orKids >= SIGODD &
                         szeKids >= MINSIGSZE & szeKids <= MAXSIGSZE)
            if (length(idx))
              x[idx]
            else
              character(0)
        })
        curCat2EntrezCond <- list()
        for (goid in names(curCat2Entrez)) {
            ## remove entrez ids that came from 
            ## SIGNIF/SIGODD/MINSIGSZE/MAXSIGSZE children
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
