\name{GeneGoCondHyperGeoTestResult-class}
\docType{class}
\alias{GeneGoCondHyperGeoTestResult-class}

\title{Class "GeneGoCondHyperGeoTestResult"} 

\description{ This class represents the results of a conditional test
  for overrepresentation of GO categories among genes in a selected gene
  set based upon the Hypergeometric distribution.  The
  \code{condHyperGeoTest} generic function returns an instance
  of the \code{GeneCategoryHyperGeoTestResult} class.

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GeneGoCondHyperGeoTestResult", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{ontology}:}{string containing the name of the GO
        ontology used for the test.  One of \code{BP}, \code{CC}, or
        \code{MF}.
      }
    \item{\code{goDag}:}{a \code{graphNEL} instance representing the
        DAG of GO terms that were tested.  Only GO terms that have at
        least one of the selected genes annotated at it are included.
      }
    \item{\code{pvalues}:}{\code{"numeric"} vector: the ordered
      p-values for each category term tested.}
    \item{\code{geneCounts}:}{\code{"integer"} vector: for each
      category term tested, the number of genes from the gene set that
      are annotated at the term.}
    \item{\code{universeCounts}:}{\code{"integer"} vector: for
      each category term tested, the number of genes from the gene
      universe that are annotated at the term.}
    \item{\code{universeMappedCount}:}{\code{"integer"} vector of length
      one: The size of the gene universe set.}
    \item{\code{geneMappedCount}:}{\code{"integer"} vector of length
      one: The size of the selected gene set.}
    \item{\code{annotation}:}{A string giving the name of the chip
      annotation data package used.}
    \item{\code{geneIds}:}{Object of class \code{"ANY"}: the input
      vector of gene identifiers passed to the testing method.}
    \item{\code{testName}:}{A string identifying the testing method
      used.}
    \item{\code{pvalue.cutoff}:}{Numeric value used a a p-value
        cutoff.  Used by the \code{show} method to count number of
        significant terms.
    }
  }
}
\section{Extends}{
Class \code{"GeneCategoryHyperGeoTestResult"}, directly.
}

\author{Seth Falcon}
\seealso{
  \code{\link[Category]{GeneCategoryHyperGeoTestResult-class}}
  \code{\link[graph]{graphNEL-class}}
  \code{\link{condHyperGeoTest}}
}
\examples{
## dummy data
res <- new("GeneGoCondHyperGeoTestResult",
           pvalues=1:10,
           geneCounts=1:10,
           universeCounts=1:10,
           universeMappedCount=100:100,
           geneMappedCount=50:50,
           annotation="hgu95av2",
                  geneIds=1:10,
           testName="CondGO",
           ontology="BP",
           goDag=new("graphNEL"))
res
}
\keyword{classes}
