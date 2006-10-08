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
