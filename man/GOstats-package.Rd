\name{GOstats-package}
\alias{GOstats-package}
\alias{GOstats}
\docType{package}
\title{
Tools for manipulating GO and microarrays.
}
\description{
A set of tools for interacting with GO and microarray
data. A variety of basic manipulation tools for graphs,
hypothesis testing and other simple calculations.
}
\details{
\tabular{ll}{
Package: \tab GOstats\cr
Version: \tab 1.7.4\cr
Date: \tab 23-08-2006\cr
biocViews: \tab Statistics, Annotation, GO, MultipleComparisons\cr
Depends: \tab graph (>= 1.9.25), GO, annotate, RBGL, xtable, Biobase,
genefilter, multtest, Category (>= 1.3.7), methods\cr
Imports: \tab methods, Category\cr
Suggests: \tab hgu95av2.db (>= 1.6.0)\cr
License: \tab Artistic\cr
}

Index:
\preformatted{
ALL                     Acute Lymphoblastic Leukemia Data from the Ritz
                        Laboratory
GOstats-defunct         Defunct Functions in GOstats Package
Ndists                  Distance matrices for the BCR/ABL and NEG
                        subgroups.
compCorrGraph           A function to compute a correlation based graph
                        from Gene Expression Data
compGdist               A function to compute the distance between
                        pairs of nodes in a graph.
dropECode               Drop GO labels for specified Evidence Codes
getEvidence             Get the Evidence codes for a set of GO terms.
getGOTerm               Functions to Access GO data.
getOntology             Get GO terms for a specified ontology
hasGOannote             Check for GO annotation
idx2dimnames            Index to Dimnames
makeGOGraph             Construct a GO Graph
notConn                 Find genes that are not connected to the
                        others.
oneGOGraph              Construct the GO graph given a set of leaves.
shortestPath            Shortest Path Analysis
simLL                   Functions to compute similarities between GO
                        graphs and also between Entrez Gene IDs based on
                        their induced GO graphs.
triadCensus             Triad Functions
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{GOstats} \tab Using GOstats (source, pdf)\cr
\code{GOusage} \tab Basic GO Usage (source, pdf)\cr
\code{GOvis} \tab Visualizing Data Using GOstats (source, pdf)\cr
}

}
\author{
R. Gentleman with contributions from S. Falcon

Maintainer: R. Gentleman <rgentlem@fhcrc.org>
}

\keyword{ package }
