\name{getOntology}
\alias{getOntology}
\title{ Get GO terms for a specified ontology }
\description{
Find the subset of GO terms, for each element of the supplied list of
bindings from a manufacturers probe.
}
\usage{
getOntology(inlist, ontology=c("MF", "BP", "CC"))
}
\arguments{
  \item{inlist}{ A list of GO associations}
  \item{ontology}{The name of the ontology you want returned. }
}
\details{
  A simple use of \code{\link{lapply}} and \code{\link{sapply}} to
  extract the terms. Only GO identifiers are returned.
}
\value{
   A list of the same length as the input list. Each element of this
   list will contain the GO term identifer for those terms that match
   the requested ontology.
}
\author{R. Gentleman}

\seealso{\code{\link{getEvidence}}, \code{\link{dropECode}} }

\examples{
 if( require("hgu95av2") ) {
     affyGO = as.list(hgu95av2GO)
     affyMF = getOntology(affyGO)
}

}
\keyword{ manip }
