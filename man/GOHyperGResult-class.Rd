\name{GOHyperGResult-class}
\docType{class}
\alias{GOHyperGResult-class}
\alias{geneCounts,GOHyperGResult-method}

\alias{goDag}
\alias{goDag,GOHyperGResult-method}

\alias{oddsRatios,GOHyperGResult-method}

\alias{expectedCounts,GOHyperGResult-method}

\alias{pvalues,GOHyperGResult-method}

\alias{geneIdUniverse,GOHyperGResult-method}

\alias{condGeneIdUniverse,GOHyperGResult-method}

\alias{htmlReport,GOHyperGResult-method}

\alias{summary,GOHyperGResult-method}

\alias{isConditional,GOHyperGResult-method}

\alias{selectedGenes}

\title{Class "GOHyperGResult"}

\description{	

  This class represents the results of a test for overrepresentation
  of GO categories among genes in a selected gene set based upon the
  Hypergeometric distribution.

  For details on extracting information from this object, be sure to
  read the accessor documentation in the Category package:
  \link[Category]{HyperGResult-accessors}.

}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GOHyperGResult", ...)}.
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
    \item{\code{pvalueCutoff}:}{Object of class \code{"numeric"}.
        The cutoff for significance used for some testing methods.
        Also used for pretty display in the \code{show} method.}

  }
}
\section{Extends}{
Class \code{"HyperGResultBase"}, directly.
}
\section{Methods}{
  \describe{
    \item{goDag}{\code{signature(r = "GOHyperGResult")}:
        return the \code{graph} instance representing the DAG of the
        GO terms that were tested.}

      \item{condGeneIdUniverse}{\code{signature(r =
        "GOHyperGResult")}: return a list named by the
      tested GO terms.  Each element of the list is a vector of the
      remaining gene identifiers (from the gene universe) annotated at
      the corresponding GO term after the identifiers from significant
      children were removed.  This is only available when the
      calculation performed used the conditional Hypergeometric
      algorithm.}

  \item{summary}{\code{signature(r = "GOHyperGResult")}: Returns a
      \code{data.frame} summarizing the test result.  Optional
      arguments \code{pvalue} and \code{categorySize} allow
      specification of maximum p-value and minimum cateogyrSize,
      respectively.  Optional argument \code{htmlLinks} is a logical
      value indicating whether to add HTML links (useful in
      conjunction with xtables print method with \code{type} set to
      \code{"html"}).}

  \item{htmlReport}{\code{signature(r = "GOHyperGResult")}: Write an
      HTML version of the table produced by the \code{summary} method.
      The path of a file to write the report to can be specified using
      the \code{file} argument.  The default is \code{file=""} which
      will cause the report to be printed to the screen.  If you wish
      to create a single report comprising multiple results you can
      set \code{append=TRUE}.  The default is \code{FALSE} (overwrite
      preexisting report file).  You can specify a string to use as an
      identifier for each table by providing a value for the
      \code{label} argument.  Additional named arguments will be
      passed to the \code{summary} method.}

  \item{description}{\code{signature(object = "GOHyperGResult")}:
      Return a string giving a one-line description of the result.}

  }


}

\author{Seth Falcon}

\seealso{
  \code{\link[Category]{GeneCategoryHyperGeoTestResult}}
}

\keyword{classes}
