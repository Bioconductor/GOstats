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
probeSetSummary(result, pvalue, categorySize, sigProbesets)
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
  \item{sigProbesets}{Optional vector of probeset IDs. See details for
      more information.}
}

\details{
  Usually the goal of doing a Fisher's exact test on a set of
  significant probesets is to find pathways or cellular activities that
  are being perturbed in an experiment. After doing the test, one
  usually gets a list of significant GO terms, and the next logical step
  might be to determine which probesets contributed to the significance
  of a certain term.

  Because the input for the Fisher's exact test consists of a vector of
  unique Entrez Gene IDs, and there may be multiple probesets that
  interrogate a particular transcript, the ouput for this function lists
  all of the probesets that map to each Entrez Gene ID, along with an
  indicator that shows which of the probesets were used as input.

  The rationale for this is that one might not be able to assume a given
  probeset actually interrogates the intended transcript, so it might be
  useful to be able to check to see what other similar probesets are
  doing.

  Because one of the first steps before running \code{hyperGTest} is to
  subset the input vectors of geneIds and universeGeneIds, any
  information about probeset IDs that interrogate the same gene
  transcript is lost. In order to recover this information, one can pass
  a vector of probeset IDs that were considered significant. This vector
  will then be used to indicate which of the probesets that map to a
  given GO term were significant in the original analysis.
}

\value{
  A \code{list} of \code{data.frame}.  Each element of the list
  corresponds to one of the GO terms (the term is provides as the name
  of the element).  Each \code{data.frame} has three columns:
  the Entrez Gene ID (\code{EntrezID}), the probe set ID
  (\code{ProbeSetID}), and a 0/1 indicator of whether the probe set ID
  was provided as part of the initial input (\code{selected})

  Note that this 0/1 indicator will only be correct if the 'geneId'
  vector used to construct the \code{GOHyperGParams} object was a named
  vector (where the names are probeset IDs), or if a vector of
  'sigProbesets' was passed to this function.
}

\author{S. Falcon and J. MacDonald}

\examples{
  ## Fake up some data
  library("hgu95av2.db")
  library("annotate")
  prbs <- ls(hgu95av2GO)[1:300]
  ## Only those with GO ids
  hasGO <- sapply(mget(prbs, hgu95av2GO), function(ids)
  if(!is.na(ids) && length(ids) > 1) TRUE else FALSE)
  prbs <- prbs[hasGO]
  prbs <- getEG(prbs, "hgu95av2")
  ## remove duplicates, but keep named vector
  prbs <- prbs[!duplicated(prbs)]
  ## do the same for universe
  univ <- ls(hgu95av2GO)[1:5000]
  hasUnivGO <- sapply(mget(univ, hgu95av2GO), function(ids)
  if(!is.na(ids) && length(ids) > 1) TRUE else FALSE)
  univ <- univ[hasUnivGO]
  univ <- unique(getEG(univ, "hgu95av2"))

  p <- new("GOHyperGParams", geneIds=prbs, universeGeneIds=univ,
  ontology="BP", annotation="hgu95av2", conditional=TRUE)
  ## this part takes time...
  if(interactive()){
    hyp <- hyperGTest(p)
    ps <- probeSetSummary(hyp, 0.05, 10)
  }
}

\keyword{manip}
\keyword{htest}
