setMethod("goDag", signature(r="GOHyperGResult"),
          function(r) r@goDag)


setMethod("pvalues", signature(r="GOHyperGResult"),
          function(r) {
              unlist(nodeData(r@goDag, attr="pvalue"))[r@pvalue.order]
              })


setMethod("oddsRatios", signature(r="GOHyperGResult"),
          function(r) {
              unlist(nodeData(r@goDag, attr="oddsRatio"))[r@pvalue.order]
              })


setMethod("expectedCounts", signature(r="GOHyperGResult"),
          function(r) {
              unlist(nodeData(r@goDag, attr="expCount"))[r@pvalue.order]
              })


entrezGeneUniverse <- function(r) {
    nodeData(r@goDag, n=nodes(r@goDag)[r@pvalue.order],
             attr="geneIds")
}


setMethod("geneIdUniverse", signature(r="GOHyperGResult"),
          function(r) {
              entrezGeneUniverse(r)
          })


setMethod("condGeneIdUniverse", signature(r="GOHyperGResult"),
          function(r) {
              if (isConditional(r))
                nodeData(r@goDag, n=nodes(r@goDag)[r@pvalue.order],
                         attr="condGeneIds")
              else
                geneIdUniverse(r)
          })


setMethod("geneCounts", signature(r="GOHyperGResult"),
          function(r) {
              sapply(condGeneIdUniverse(r), function(x) {
                  sum(r@geneIds %in% x)
              })
          })


## setMethod("condGeneCounts", signature(r="GOHyperGResult"),
##           function(r) {
##               sapply(condGeneIdUniverse(r), function(x) {
##                   sum(r@geneIds %in% x)
##               })
##           })


setMethod("universeCounts", signature(r="GOHyperGResult"),
          function(r) {
              sapply(entrezGeneUniverse(r), length)
          })


setMethod("universeMappedCount", signature(r="GOHyperGResult"),
          function(r) {
              length(unique(unlist(entrezGeneUniverse(r))))
          })


setMethod("isConditional", signature(r="GOHyperGResult"),
          function(r) r@conditional)


setMethod("description", signature(object="GOHyperGResult"),
          function(object) {
              cond <- "Conditional"
              if (!isConditional(object))
                cond <- ""
              desc <- paste("Gene to %s", cond, "test for %s-representation")
              desc <- sprintf(desc,
                              paste(testName(object), collapse=" "),
                              testDirection(object))
              desc
          })


selectedGenes <- function(r, id) {
    ## FIXME: make me a method!
    ans <- geneIdUniverse(r)[id]
    ans <- lapply(ans, intersect, geneIds(r))
    ans
}


sigCategories <- function(res, p) {
    ## FIXME: make me a method!
    if (missing(p))
      p <- pvalueCutoff(res)
    pv <- pvalues(res)
    goIds <- names(pv[pv < p])
    goIds
}


inducedTermGraph <- function(r, id, children=TRUE, parents=TRUE,
                             ...) {
    if (!children && !parents)
      stop("children and parents can't both be FALSE")
    ## XXX: should use more structure here
    goName <- paste(testName(r), collapse="")
    goKidsEnv <- get(paste(goName, "CHILDREN", sep=""))
    goParentsEnv <- get(paste(goName, "PARENTS", sep=""))
    goIds <- character(0)

    wantedNodes <- id
    ## children
    if (children) {
        wantedNodes <- c(wantedNodes,
                         unlist(edges(goDag(r))[id], use.names=FALSE))
    }
    ## parents
    g <- reverseEdgeDirections(goDag(r))
    if (parents) {
        wantedNodes <- c(wantedNodes,
                         unlist(edges(g)[id], use.names=FALSE))
    }
    wantedNodes <- unique(wantedNodes)
    g <- subGraph(wantedNodes, g)

    ## expand
    if (children) {
        for (goid in id) {
            kids <- unique(goKidsEnv[[goid]])
            for (k in kids) {
                if (!(k %in% nodes(g))) {
                    g <- addNode(k, g)
                    g <- addEdge(k, goid, g)
                }
            }
        }
    }
    if (parents) {
        for (goid in id) {
            elders <- unique(goKidsEnv[[goid]])
            for (p in elders) {
                if (!(p %in% nodes(g))) {
                    g <- addNode(p, g)
                    g <- addEdge(goid, p, g)
                }
            }
        }
    }
    g
}


plotGOTermGraph <- function(g, r=NULL, add.counts=TRUE,
                            max.nchar=20,...) {
    n <- nodes(g)
    termLab <- substr(sapply(mget(n, GOTERM), Term), 0, max.nchar)
    ncolors <- rep("red", length(n))
    if (!is.null(r) && add.counts) {
        resultTerms <- names(pvalues(r))
        ncolors <- ifelse(n %in% resultTerms, "lightgray", "white")
        counts <- sapply(n, function(x) {
            if (x %in% resultTerms) {
                paste(geneCounts(r)[x], "/",
                      universeCounts(r)[x],
                      sep="")
            } else {
                ""
            }
        })
        nlab <- paste(termLab, counts)
    } else {
        nlab <- termLab
    }
    nattr <- makeNodeAttrs(g,
                           label=nlab,
                           shape="ellipse",
                           fillcolor=ncolors,
                           fixedsize=FALSE)
    plot(g, ..., nodeAttrs=nattr)
}


termGraphs <- function(r, id=NULL, pvalue=NULL, use.terms=TRUE) {
    if (!is.null(id) && !is.null(pvalue))
      warning("ignoring pvalue arg since GO IDs where specified")
    if (missing(pvalue) || is.null(pvalue))
      pvalue <- pvalueCutoff(r)
    if (is.null(id))
      goids <- sigCategories(r, pvalue)
    else
      goids <- id
    g <-  reverseEdgeDirections(subGraph(goids, goDag(r)))
    cc <- connectedComp(g)
    sapply(cc, subGraph, g)
}
