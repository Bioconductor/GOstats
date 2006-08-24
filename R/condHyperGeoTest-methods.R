removeLengthZero <- function(x) {
    wanted <- sapply(x, function(z) length(z) > 0)
    x[wanted]
}

setMethod("condHyperGeoTest",
          signature(p="GeneCategoryHyperGeoTestParams"), 
          function(p) {
              origGeneIds <- p@geneIds
              p@universeGeneIds <- universeBuilder(p)
              selected <- intersect(p@geneIds, p@universeGeneIds)
              p@geneIds <- selected
              cat2Entrez <- categoryToEntrezBuilder(p)
              ## build the GO graph for the relevant GO nodes
              goIds <- names(cat2Entrez)
              cat("found", length(goIds), "GO ids to test\n")
              ## first, reduce GOBPCHILDREN to our current set
              ## FIXME: need a method here
              gobpkids <- eapply(GOBPCHILDREN, function(x) {
                  intersect(x, goIds)
              })
              gobpkids <- l2e(gobpkids)
              goDag <- GOGraph(goIds, gobpkids)
              nodeDataDefaults(goDag, "pvalue") <- 1
              ## now iterate leaves first doing tests and conditioning
              ## on all significant children.
              needsProc <- goDag
              complete <- character(0)
              SIGNIF <- p@pvalue.cutoff
              while (length(nodes(needsProc))) {
                  numKids <- sapply(edges(needsProc), length)
                  noKids <- names(numKids[numKids == 0])
                  curCat2Entrez <- cat2Entrez[noKids]
                  curCatKids <- edges(goDag)[names(curCat2Entrez)]
                  curCatKids <- removeLengthZero(curCatKids)
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
                  numFound <- sapply(curCat2Entrez, 
                                     function(x) sum(selected %in% x))
                  numDrawn <- length(selected)
                  ## num white in urn
                  numAtCat <- sapply(curCat2Entrez, length)
                  ## num black in urn
                  numNotAtCat <- length(p@universeGeneIds) - numAtCat
                  ## take the -1 because we want evidence for as extreme or more.
                  pvals <- phyper(numFound-1, numAtCat, numNotAtCat, numDrawn,
                                  lower.tail=FALSE)
                  ## store the pvals, mark these nodes as complete,
                  ## then compute the next set of nodes to do.  
                  noKids <- names(curCat2Entrez)
                  ## drop names on pvals to avoid weird names upon unlisting
                  nodeData(goDag, n=noKids, attr="pvalue") <- as.numeric(pvals)
                  complete <- c(complete, noKids)
                  hasKids <- names(numKids[numKids > 0])
                  needsProc <- subGraph(hasKids, needsProc)
              }
              pvals <- unlist(nodeData(goDag, attr="pvalue"))
              ord <- order(pvals)
              new("GeneGoCondHyperGeoTestResult",
                  pvalues=pvals[ord],
                  geneCounts=numFound[ord],
                  universeCounts=numAtCat[ord],
                  universeMappedCount=length(p@universeGeneIds),
                  geneMappedCount=length(p@geneIds),
                  annotation=p@annotation,
                  geneIds=origGeneIds,
                  testName=p@categoryName,
                  ontology=p@ontology,
                  goDag=goDag,
                  pvalue.cutoff=p@pvalue.cutoff)
          })
