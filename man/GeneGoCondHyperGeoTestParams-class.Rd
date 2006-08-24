\name{GeneGoCondHyperGeoTestParams-class}
\docType{class}
\alias{GeneGoCondHyperGeoTestParams-class}

\title{Class "GeneGoCondHyperGeoTestParams"}
\description{	
  A parameter class for representing all parameters needed for running
  the \code{condHyperGeoTest} method with one of the GO
  ontologies (BP, CC, MF) as the category.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GeneGoCondHyperGeoTestParams", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{pvalue.cutoff}:}{\code{"numeric"} value that will be
        used to determine significance during calculations.}
    \item{\code{ontology}:}{A string specifying the GO ontology to use.
      Must be one of "BP", "CC", or "MF".}
    \item{\code{geneIds}:}{Object of class \code{"ANY"}: A vector of
      gene identifiers.  Numeric and character vectors are probably the
      only things that make sense.  These are the gene ids for the
      selected gene set.}
    \item{\code{universeGeneIds}:}{Object of class \code{"ANY"}: A
      vector of gene ids in the same format as \code{geneIds} defining a
      subset of the gene ids on the chip that will be used as the
      universe for the hypergeometric calculation.  If this is
      \code{NULL} or has length zero, then all gene ids on the chip will
      be used.}
    \item{\code{annotation}:}{A string giving the name of the
      annotation data package for the chip used to generate the data.}
    \item{\code{cateogrySubsetIds}:}{Object of class \code{"ANY"}:
      If the test method supports it, can be used to specify a subset of
      category ids to include in the test instead of all possible
      category ids.}
    \item{\code{categoryName}:}{A string describing the category.
      Usually set automatically by subclasses.  For example "GO".}           
  }
}
\section{Extends}{
Class \code{"GeneGoHyperGeoTestParams"}, directly.
Class \code{"GeneCategoryHyperGeoTestParams"}, by class \code{"GeneGoHyperGeoTestParams"}.
}
\section{Methods}{
No methods defined with class "GeneGoCondHyperGeoTestParams" in the signature.
}

\author{Seth Falcon}

}
\examples{
## Extract a random subset of Entrez Gene IDs as
## the gene universe and another random sample of 
## "selected" genes.
    set.seed(344)
    probeIds <- ls(hgu95av2LOCUSID)
    randProbeIds <- sample(probeIds, 500)
    entrezUniverse <- mget(randProbeIds, hgu95av2LOCUSID, ifnotfound=NA)
    entrezUniverse <- entrezUniverse[!is.na(entrezUniverse)]
    selectedEntrezIds <- sample(entrezUniverse, 30)

## Create a new parameter instance
    params <- new("GeneGoCondHyperGeoTestParams",
                  geneIds=selectedEntrezIds, 
                  universeGeneIds=entrezUniverse,
                  annotation="hgu95av2", 
                  ontology="BP",
                  pvalue.cutoff=0.01)
    params

}
\keyword{classes}
