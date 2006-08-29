\name{GeneGoHyperGeoTestResult-class}
\docType{class}
\alias{GeneGoHyperGeoTestResult-class}
\alias{geneCounts,GeneGoHyperGeoTestResult-method}

\alias{goDag}
\alias{goDag,GeneGoHyperGeoTestResult-method}

\alias{pvalues,GeneGoHyperGeoTestResult-method}
\alias{universeCounts,GeneGoHyperGeoTestResult-method}
\alias{universeMappedCount,GeneGoHyperGeoTestResult-method}

\title{Class "GeneGoHyperGeoTestResult"}

\description{	

  This class represents the results of a test for overrepresentation
  of GO categories among genes in a selected gene set based upon the
  Hypergeometric distribution.

}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GeneGoHyperGeoTestResult", ...)}.
}

\section{Slots}{
	 \describe{
    \item{\code{goDag}:}{Object of class \code{"graph"} representing
        the DAG of GO terms tested.}
    \item{\code{pvalue.order}:}{Object of class \code{"integer"}.
        The sort order of the computed p-values.}
    \item{\code{annotation}:}{Object of class \code{"character"}.  The
      name of the annotation data package used in the analysis.}
    \item{\code{geneIds}:}{Object of class \code{"ANY"}.  The
        intersection of the gene identifiers given as input and the
        computed gene universe.}
    \item{\code{testName}:}{Object of class \code{"character"}.
        Identifies the testing method used to produce this result instance.}
    \item{\code{pvalue.cutoff}:}{Object of class \code{"numeric"}.
        The cutoff for significance used for some testing methods.
        Also used for pretty display in the \code{show} method.}

  }
}
\section{Extends}{
Class \code{"GeneCategoryHyperGeoTestResultBase"}, directly.
}
\section{Methods}{
  \describe{
  \item{geneCounts}{\code{signature(r =
        "GeneGoHyperGeoTestResult")}: return an \code{"integer"}
      vector: for each category term tested, the number of genes from
      the gene set that are annotated at the term.}

    \item{goDag}{\code{signature(r = "GeneGoHyperGeoTestResult")}:
        return the \code{graph} instance representing the DAG of the
        GO terms that were tested.}

 \item{pvalues}{\code{signature(r =
        "GeneGoHyperGeoTestResult")}: return a \code{"numeric"}
      vector: the ordered p-values for each category term tested.}

  \item{universeCounts}{\code{signature(r =
        "GeneGoHyperGeoTestResult")}: return an \code{"integer"}
      vector: for each category term tested, the number of genes from
      the gene universe that are annotated at the term.}

  \item{universeMappedCount}{\code{signature(r =
        "GeneGoHyperGeoTestResult")}: return an \code{"integer"}
      vector of length one giving the size of the gene universe set. }

  }
}

\author{Seth Falcon}

\seealso{
  \code{\link[Category]{GeneCategoryHyperGeoTestResult}}
}

\keyword{classes}
