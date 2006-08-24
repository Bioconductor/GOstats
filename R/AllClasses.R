setClass("GeneGoCondHyperGeoTestParams",
         representation(pvalue.cutoff="numeric"),
         contains="GeneGoHyperGeoTestParams",
         prototype=prototype(pvalue.cutoff=0.01))


setClass("GeneGoCondHyperGeoTestResult", 
         contains="GeneCategoryHyperGeoTestResult",
         representation=representation(
           ontology="character",
           goDag="graphNEL"))
