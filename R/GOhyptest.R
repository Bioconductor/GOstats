resultToGOHyperG <- function(res, origIds) {
    go2Affy <- mget(names(pvalues(res)),
                    Category:::getDataEnv("GO2ALLPROBES", annotation(res)))
    list(pvalues=pvalues(res),
         goCounts=universeCounts(res),
         intCounts=geneCounts(res),
         numLL=universeMappedCount(res),
         numInt=geneMappedCount(res),
         chip=annotation(res),
         intLLs=origIds,
         go2Affy=go2Affy)
}


GOHyperG <- function(x, lib, what="MF", universe=NULL)
{
    if (missing(universe) || is.null(universe))
      universe <- character(0)
    params <- new("GeneGoHyperGeoTestParams",
                  geneIds=x,
                  universeGeneIds=universe,
                  annotation=lib,
                  ontology=what)
    if (missing(lib) || !is.character(lib))
      stop("argument ", sQuote("lib"), " must be character")

    res <- geneCategoryHyperGeoTest(params)
    resultToGOHyperG(res, origIds=x)
}


GOKEGGHyperG <- function (x, lib = "hgu95av2", what = "MF", universe=NULL)
{
    .Defunct("geneKeggHyperGeoTest", "Category")
}
