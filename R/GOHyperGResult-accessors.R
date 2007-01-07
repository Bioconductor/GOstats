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


plotTermGraphs <- function(r, pvalue=NULL, use.terms=TRUE,
                           dev="x11", ...) {
    if (missing(pvalue) || is.null(pvalue))
      pvalue <- pvalueCutoff(r)
    goids <- sigCategories(r, pvalue)
    g <-  reverseEdgeDirections(subGraph(goids, goDag(r)))
    cc <- connectedComp(g)
    if (!interactive())
      dev <- getOption("device")
    openMany <- TRUE
    if (!(dev %in% c("x11", "X11", "quartz", "windows", "win.graph")))
      openMany <- FALSE
    theDev <- get(dev)
    if (!openMany)
      theDev(...)
    for (n in cc) {
        sg <- subGraph(n, g)
        if (use.terms) {
            termLab <- sapply(mget(nodes(sg), GOTERM), Term)
        } else {
            termLab <- sub("^GO:", "", nodes(sg))
        }
        nlab <- paste(termLab, " (",
                      geneCounts(r)[n], "/",
                      universeCounts(r)[n],
                      ")", sep="")
        if (length(n) == 1) {
            message("Skipping singleton component: ",
                    nlab)
            next
        }
        nattr <- makeNodeAttrs(sg,
                               label=nlab,
                               shape="ellipse",
                               fixedsize=FALSE)
        if (openMany) 
          theDev(...)
        plot(sg, nodeAttrs=nattr)
    }
    if (!openMany)
      dev.off()
}
