\name{OBOHyperGResult-class}
\docType{class}
\alias{OBOHyperGResult-class}
\alias{geneCounts,OBOHyperGResult-method}

\alias{goDag,OBOHyperGResult-method}

\alias{oddsRatios,OBOHyperGResult-method}

\alias{expectedCounts,OBOHyperGResult-method}

\alias{pvalues,OBOHyperGResult-method}

\alias{geneIdUniverse,OBOHyperGResult-method}

\alias{htmlReport,OBOHyperGResult-method}

\alias{summary,OBOHyperGResult-method}

\alias{conditional,OBOHyperGResult-method}

\title{Class "OBOHyperGResult"}

\description{	

  This class represents the results of a test for overrepresentation
  of OBO categories among genes in a selected gene set based upon the
  Hypergeometric distribution.

  For details on extracting information from this object, be sure to
  read the accessor documentation in the Category package:
  \link[Category]{HyperGResult-accessors}.

}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("OBOHyperGResult", ...)}.
}

\section{Slots}{
	 \describe{
    \item{\code{goDag}:}{Object of class \code{"graph"} representing
        the DAG of OBO terms tested.}
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
    \item{\code{conditional}:}{A logical indicating whether the
      calculation should condition on the OBO structure.}
    \item{\code{testDirection}:}{A string which can be either "over"
      or "under". This determines whether the test performed detects
      over or under represented OBO terms.}
  }
}

\section{Extends}{
Class \code{"HyperGResultBase"}, directly.
}
\section{Methods}{
  \describe{
    \item{goDag}{\code{signature(r = "OBOHyperGResult")}:
        return the \code{graph} instance representing the DAG of the
        OBO terms that were tested.}


  \item{summary}{\code{signature(r = "OBOHyperGResult")}: Returns a
      \code{data.frame} summarizing the test result.  Optional
      arguments \code{pvalue} and \code{categorySize} allow
      specification of maximum p-value and minimum categorySize,
      respectively.  Optional argument \code{htmlLinks} is a logical
      value indicating whether to add HTML links (useful in
      conjunction with xtables print method with \code{type} set to
      \code{"html"}).}

  \item{htmlReport}{\code{signature(r = "OBOHyperGResult")}: Write an
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

  \item{description}{\code{signature(object = "OBOHyperGResult")}:
      Return a string giving a one-line description of the result.}

  }


}

\author{Robert Castelo}

\seealso{
  \link[Category]{HyperGResult-accessors}
}

\keyword{classes}
