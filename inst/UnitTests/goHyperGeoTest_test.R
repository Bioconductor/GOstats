makeSimpleGeneGoHyperGeoTestParams <- function() {
    set.seed(344)
    probeIds <- ls(hgu95av2LOCUSID)
    randProbeIds <- sample(probeIds, 500)
##     entrezUniverse <- unlist(mget(randProbeIds, hgu95av2LOCUSID,
##                                   ifnotfound=NA))
    ## This is "wrong", should unlist, but the code
    ## should catch/correct it.  The right way is above.
    entrezUniverse <- mget(randProbeIds, hgu95av2LOCUSID,
                           ifnotfound=NA)

    entrezUniverse <- entrezUniverse[!is.na(entrezUniverse)]
    selectedEntrezIds <- sample(entrezUniverse, 30)
    params <- new("GeneGoHyperGeoTestParams",
                  geneIds=selectedEntrezIds, 
                  universeGeneIds=entrezUniverse,
                  annotation="hgu95av2", 
                  ontology="BP",
                  pvalue.cutoff=0.05,
                  conditional=TRUE,
                  test.direction="over")
    params
}
    

test_goHyperGeoTest_regression1 <- function() {
    p <- makeSimpleGeneGoHyperGeoTestParams()
    res <- goHyperGeoTest(p)
    checkEquals(9, sum(pvalues(res) < res@pvalue.cutoff))
    first5 <- c("GO:0009057", "GO:0019320", "GO:0006096",
                "GO:0006006", "GO:0046164")
    checkEquals(first5, names(pvalues(res)[1:5]))
    checkEquals(166, length(universeCounts(res)))
    checkEquals(253, numEdges(goDag(res)))
}


test_goHyperGeoTest_regression2 <- function() {
    p <- makeSimpleGeneGoHyperGeoTestParams()
    p@conditional <- FALSE
    res <- goHyperGeoTest(p)
    checkEquals(18, sum(pvalues(res) < res@pvalue.cutoff))
    
    res2 <- geneCategoryHyperGeoTest(p)
    checkEquals(pvalues(res), pvalues(res2))
    checkEquals(geneCounts(res), geneCounts(res2))
    checkEquals(universeCounts(res), universeCounts(res2))
    checkEquals(universeMappedCount(res), universeMappedCount(res2))
    checkEquals(geneMappedCount(res), geneMappedCount(res2))
    checkEquals(annotation(res), annotation(res2))
    checkEquals(testName(res), testName(res2))
}
