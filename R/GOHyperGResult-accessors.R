setMethod("goDag", signature(r="GOHyperGResult"),
          function(r) r@goDag)

orderedAttr <- function(r, attr) {
    unlist(nodeData(r@goDag, attr=attr))[r@pvalue.order]
}

setMethod("pvalues", signature(r="GOHyperGResult"),
          function(r) orderedAttr(r, "pvalue"))

setMethod("oddsRatios", signature(r="GOHyperGResult"),
          function(r) orderedAttr(r, "oddsRatio"))

setMethod("expectedCounts", signature(r="GOHyperGResult"),
          function(r) orderedAttr(r, "expCount"))


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

setMethod("isConditional", signature(r="GOHyperGResult"),
          function(r) r@conditional)


selectedGenes <- function(r, id=NULL) {
    .Deprecated("geneIdsByCategory", package="Category")
    geneIdsByCategory(r, id)
}


inducedTermGraph <- function(r, id, children=TRUE, parents=TRUE) {
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

    ## expand; add children and/or parents that are not present in g,
    ## but are definedin the GO data.
    if (children) {
        for (goid in id) {
            kids <- unique(goKidsEnv[[goid]])
            for (k in kids) {
                if (is.na(k)) next
                if (!(k %in% nodes(g))) {
                    g <- addNode(k, g)
                    g <- addEdge(k, goid, g)
                }
            }
        }
    }
    if (parents) {
        for (goid in id) {
            elders <- unique(goParentsEnv[[goid]])
            for (p in elders) {
                if (is.na(p)) next
                if (!(p %in% nodes(g))) {
                    g <- addNode(p, g)
                    g <- addEdge(goid, p, g)
                }
            }
        }
    }
    g
}


## FIXME: perhpas it doesn't make sense to exclude the untestable GO terms.
## maybe it would be better to keep them in as it will be less confusing?
plotGOTermGraph <- function(g, r=NULL, add.counts=TRUE,
                            max.nchar=20,
                            node.colors=c(sig="lightgray", not="white"),
                            ...) {
    if (!require("Rgraphviz", quietly=TRUE))
      stop("The Rgraphviz package is required for this feature")
    termLab <- n <- nodes(g)
    if (!is.null(max.nchar))
      termLab <- sapply(termLab, substr, 1L, max.nchar, USE.NAMES=FALSE)
    ncolors <- rep(node.colors["not"], length(n))
    if (!is.null(r) && add.counts) {
        if (is.null(names(node.colors)) ||
            !all(c("sig", "not") %in% names(node.colors)))
          stop(paste("invalid node.colors arg:",
                     "must have named elements 'sig' and 'not'"))
        resultTerms <- names(pvalues(r))
        ncolors <- ifelse(n %in% sigCategories(r), node.colors["sig"],
                          node.colors["not"])
        counts <- sapply(n, function(x) {
            if (x %in% resultTerms) {
                paste(geneCounts(r)[x], "/",
                      universeCounts(r)[x],
                      sep="")
            } else {
                "0/??"
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
    g <- subGraph(goids, goDag(r))
    if (use.terms)
      nodes(g) <- as.character(sapply(mget(nodes(g), GOTERM), Term))
    g <-  reverseEdgeDirections(g)
    cc <- connectedComp(g)
    sapply(cc, subGraph, g)
}
