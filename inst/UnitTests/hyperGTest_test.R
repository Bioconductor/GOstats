library("hgu95av2.db")

makeSimpleGOHyperGParams <- function() {
    set.seed(344)
    probeIds <- ls(hgu95av2ENTREZID)
    randProbeIds <- sample(probeIds, 500)
    ## This is "wrong", should unlist, but the code
    ## should catch/correct it.  The right way is to
    ## unlist the mget result.
    entrezUniverse <- mget(randProbeIds, hgu95av2ENTREZID,
                           ifnotfound=NA)
    entrezUniverse <- entrezUniverse[!is.na(entrezUniverse)]
    selectedEntrezIds <- sample(entrezUniverse, 30)
    params <- new("GOHyperGParams",
                  geneIds=selectedEntrezIds, 
                  universeGeneIds=entrezUniverse,
                  annotation="hgu95av2", 
                  ontology="BP",
                  pvalueCutoff=0.05,
                  conditional=TRUE,
                  testDirection="over")
    params
}
    


test_hyperGTest_regression2 <- function() {
    p <- makeSimpleGOHyperGParams()
    p@conditional <- FALSE
    res <- hyperGTest(p)

    ## Verify result is same using geneCategoryHyperGeoTest
    ## from Category
    res2 <- Category::hyperGTest(p)
    checkEquals(pvalues(res), pvalues(res2))
    checkEquals(geneCounts(res), geneCounts(res2))
    checkEquals(universeCounts(res), universeCounts(res2))
    checkEquals(universeMappedCount(res), universeMappedCount(res2))
    checkEquals(geneMappedCount(res), geneMappedCount(res2))
    checkEquals(annotation(res), annotation(res2))
    checkEquals(testName(res), testName(res2))

    ## Regression tests
    checkEquals("hgu95av2", annotation(res))
    checkEquals(c("GO", "BP"), testName(res))
}
