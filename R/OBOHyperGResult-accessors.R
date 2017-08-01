## this is for OBO
## the function orderedAttr() is defined
## in GOHyperGResults-accessors.R

setMethod("pvalues", signature(r="OBOHyperGResult"),
          function(r) orderedAttr(r, "pvalue"))

setMethod("oddsRatios", signature(r="OBOHyperGResult"),
          function(r) orderedAttr(r, "oddsRatio"))

setMethod("expectedCounts", signature(r="OBOHyperGResult"),
          function(r) orderedAttr(r, "expCount"))

## the function entrezGeneIdUniverse() is defined
## in GOHyperGResults-accessors.R

setMethod("geneIdUniverse", signature(r="OBOHyperGResult"),
          function(r, cond=TRUE) {
            if (cond && conditional(r))
              nodeData(r@goDag, n=nodes(r@goDag)[r@pvalue.order],
                       attr="condGeneIds")
            else
              entrezGeneUniverse(r)
          })

setMethod("summary", signature(object="OBOHyperGResult"),
          function(object, pvalue=pvalueCutoff(object),
                   categorySize=NULL) {
            df <- callNextMethod(object=object, pvalue=pvalue,
                                 categorySize=categorySize)
            if (nrow(df) == 0) {
              df$Term <- character(0)
              return(df)
            }
            oboIds <- df[[1]]
            df$Term <- object@gscDescriptions[oboIds]
            df
          })
