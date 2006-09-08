setClass("GeneGoHyperGeoTestResult",
         contains="GeneCategoryHyperGeoTestResultBase",
         representation=representation(
           goDag="graph",
           pvalue.order="integer"),
         prototype=prototype(
           testName="GO",
           pvalueCutoff=0.01,
           goDag=new("graphNEL")))


