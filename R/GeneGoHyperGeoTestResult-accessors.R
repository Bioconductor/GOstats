setMethod("goDag", signature(r="GeneGoHyperGeoTestResult"),
          function(r) r@goDag)


setMethod("pvalues", signature(r="GeneGoHyperGeoTestResult"),
          function(r) {
              unlist(nodeData(r@goDag, attr="pvalue"))[r@pvalue.order]
              })


setMethod("oddsRatios", signature(r="GeneGoHyperGeoTestResult"),
          function(r) {
              unlist(nodeData(r@goDag, attr="oddsRatio"))[r@pvalue.order]
              })


setMethod("expectedCounts", signature(r="GeneGoHyperGeoTestResult"),
          function(r) {
              unlist(nodeData(r@goDag, attr="expCount"))[r@pvalue.order]
              })


entrezGeneUniverse <- function(r) {
    nodeData(r@goDag, n=nodes(r@goDag)[r@pvalue.order],
             attr="geneIds")
}


setMethod("geneIdUniverse", signature(r="GeneGoHyperGeoTestResult"),
          function(r) {
              entrezGeneUniverse(r)
          })


setMethod("condGeneIdUniverse", signature(r="GeneGoHyperGeoTestResult"),
          function(r) {
              if ("condGeneIds" %in% names(nodeDataDefaults(r@goDag)))
                nodeData(r@goDag, n=nodes(r@goDag)[r@pvalue.order],
                         attr="condGeneIds")
              else
                NA
          })


setMethod("geneCounts", signature(r="GeneGoHyperGeoTestResult"),
          function(r) {
              sapply(entrezGeneUniverse(r), function(x) {
                  sum(r@geneIds %in% x)
              })
          })


setMethod("universeCounts", signature(r="GeneGoHyperGeoTestResult"),
          function(r) {
              sapply(entrezGeneUniverse(r), length)
          })


setMethod("universeMappedCount", signature(r="GeneGoHyperGeoTestResult"),
          function(r) {
              length(unique(unlist(entrezGeneUniverse(r))))
          })
