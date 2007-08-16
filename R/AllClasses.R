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

setClass("KEGGHyperGResult",
         contains="HyperGResult")

setClass("PFAMHyperGResult",
         contains="HyperGResult")


