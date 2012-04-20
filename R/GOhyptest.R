resultToGOHyperG <- function(res, origIds) {
    go2Affy <- mget(names(pvalues(res)),
                    getAnnMap("GO2ALLPROBES", annotation(res)))
    list(pvalues=pvalues(res),
         goCounts=universeCounts(res),
         intCounts=geneCounts(res),
         numLL=universeMappedCount(res),
         numInt=geneMappedCount(res),
         chip=annotation(res),
         intLLs=origIds,
         go2Affy=go2Affy)
}
