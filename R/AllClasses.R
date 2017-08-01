setClass("GOHyperGResult",
         contains="HyperGResultBase",
         representation=representation(
           goDag="graph",
           pvalue.order="integer",
           conditional="logical"),
         prototype=prototype(
           testName="GO",
           pvalueCutoff=0.01,
           goDag=new("graphNEL")))

## this is for OBO
setClass("OBOHyperGResult",
         contains="HyperGResultBase",
         representation=representation(
           goDag="graph",
           gscDescriptions="character",
           pvalue.order="integer",
           conditional="logical"),
         prototype=prototype(
           testName="GO",
           pvalueCutoff=0.01,
           goDag=new("graphNEL")))
