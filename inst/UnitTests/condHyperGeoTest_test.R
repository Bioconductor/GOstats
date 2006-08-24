require("hgu95av2") || stop("hgu95av2 is needed for unit tests")

makeSimpleGeneGoHyperGeoTestParams <- function() {
    set.seed(344)
    probeIds <- ls(hgu95av2LOCUSID)
    randProbeIds <- sample(probeIds, 500)
    entrezUniverse <- mget(randProbeIds, hgu95av2LOCUSID, ifnotfound=NA)
    entrezUniverse <- entrezUniverse[!is.na(entrezUniverse)]
    selectedEntrezIds <- sample(entrezUniverse, 30)
    params <- new("GeneGoCondHyperGeoTestParams",
                  geneIds=selectedEntrezIds, 
                  universeGeneIds=entrezUniverse,
                  annotation="hgu95av2", 
                  ontology="BP",
                  pvalue.cutoff=0.05)
    params
}
    

test_condHyperGeoTest_regression <- function() {
    p <- makeSimpleGeneGoHyperGeoTestParams()
    res <- condHyperGeoTest(p)
    checkEquals(18, sum(res@pvalues < res@pvalue.cutoff))
    first5 <- c("GO:0043170", "GO:0044265", "GO:0009057",
                "GO:0044248", "GO:0019320")
    checkEquals(first5, names(res@pvalues[1:5]))
    checkEquals(166, numNodes(res@goDag))
    checkEquals(253, numEdges(res@goDag))
}
