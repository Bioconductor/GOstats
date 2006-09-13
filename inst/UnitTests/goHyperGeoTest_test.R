makeSimpleGOHyperGParams <- function() {
    set.seed(344)
    probeIds <- ls(hgu95av2LOCUSID)
    randProbeIds <- sample(probeIds, 500)
    ## This is "wrong", should unlist, but the code
    ## should catch/correct it.  The right way is to
    ## unlist the mget result.
    entrezUniverse <- mget(randProbeIds, hgu95av2LOCUSID,
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
    

test_goHyperGeoTest_regression1 <- function() {
    p <- makeSimpleGOHyperGParams()
    res <- goHyperGeoTest(p)
    checkEquals(9, sum(pvalues(res) < res@pvalueCutoff))
    first5 <- c("GO:0009057", "GO:0019320", "GO:0006096",
                "GO:0006006", "GO:0046164")
    checkEquals(first5, names(pvalues(res)[1:5]))
    checkEquals(166, length(universeCounts(res)))
    checkEquals(253, numEdges(goDag(res)))
}


test_goHyperGeoTest_regression2 <- function() {
    p <- makeSimpleGOHyperGParams()
    p@conditional <- FALSE
    res <- goHyperGeoTest(p)
    checkEquals(18, sum(pvalues(res) < res@pvalueCutoff))

    ## Verify result is same using geneCategoryHyperGeoTest
    ## from Category
    res2 <- geneCategoryHyperGeoTest(p)
    checkEquals(pvalues(res), pvalues(res2))
    checkEquals(geneCounts(res), geneCounts(res2))
    checkEquals(universeCounts(res), universeCounts(res2))
    checkEquals(universeMappedCount(res), universeMappedCount(res2))
    checkEquals(geneMappedCount(res), geneMappedCount(res2))
    checkEquals(annotation(res), annotation(res2))
    checkEquals(testName(res), testName(res2))

    ## Regression tests
        pvals <- round(c(0.01596240, 0.01825930, 0.02194761), 3)
    names(pvals) <- c("GO:0043170", "GO:0044265", "GO:0009057") 
    checkEquals(pvals, round(pvalues(res)[1:3], 3))
    
    gcounts <- c(13, 4, 4)
    names(gcounts) <- c("GO:0043170", "GO:0044265", "GO:0009057")
    checkEquals(gcounts, geneCounts(res)[1:3])

    ucounts <- c(134, 19, 20)
    names(ucounts) <- c("GO:0043170", "GO:0044265", "GO:0009057")
    checkEquals(ucounts, universeCounts(res)[1:3])
    
    checkEquals(381, universeMappedCount(res))
    checkEquals(22, geneMappedCount(res))
    checkEquals("hgu95av2", annotation(res))
    checkEquals(c("GO", "BP"), testName(res))

    ## Test conversion to legacy result

    ## Check conversion to GOHyperG output format
    ## We use a serialized result and test only the first
    ## 10 entries (otherwise the saved result is too large).
    ghg <- GOstats:::resultToGOHyperG(res, p@geneIds)
    load(system.file("UnitTests/ghgans10.rda",
                     package="GOstats"))
    ghg10_expected <- ghgans10 ## from load
    ghg10 <- lapply(ghg, function(x) if (length(x) > 11) x[1:10] else x)
    checkEquals(ghg10_expected, ghg10)

}
