\name{probeSetSummary}
\alias{probeSetSummary}

\title{Summarize Probe Sets Associated with a hyperGTest Result}
\description{
  Given the result of a \code{hyperGTest} run (an instance of
  \code{GOHyperGResult}), this function lists all Probe Set IDs
  associated with the selected Entrez IDs annotated at each
  significant GO term in the test result.
}
\usage{
probeSetSummary(result, pvalue, categorySize)
}

\arguments{
  \item{result}{A \code{GOHyperGResult} instance.  This is the output   
      of the \code{hyperGTest} function when testing the GO category.} 

  \item{pvalue}{Optional p-value cutoff.  Only results for GO terms
      with a p-value less than the specified value will be returned.
      If omitted, \code{pvalueCutoff(result)} is used.}
  \item{categorySize}{Optional minimum size (number of annotations)
      for the GO terms.  Only results for GO terms with
      \code{categorySize} or more annotations will be returned.  If
      omitted, no category size criteria will be used.}
}

\value{
  A \code{list} of \code{data.frame}.  Each element of the list
  corresponds to one of the GO terms (the term is provides as the name
  of the element).  Each \code{data.frame} has three columns:
  the Entrez Gene ID (\code{EntrezID}), the probe set ID
  (\code{ProbeSetID}), and a 0/1 indicator of whether the probe set ID
  was provided as part of the initial input (\code{selected}).
}

\author{S. Falcon and J. MacDonald}

\keyword{manip}
\keyword{htest}
