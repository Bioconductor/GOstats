setClass("GeneGoHyperGeoTestResult",
         contains="GeneCategoryHyperGeoTestResultBase",
         representation=representation(
           goDag="graph",
           pvalue.order="integer"),
         prototype=prototype(
           testName="GO",
           pvalue.cutoff=0.01,
           goDag=new("graphNEL")))


