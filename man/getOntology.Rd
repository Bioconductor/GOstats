\name{getOntology}
\alias{getOntology}
\title{ Get GO terms for a specified ontology }
\description{
Find the subset of GO terms for the specified ontology, for each element 
of the supplied list of associations. The input list is typically from
one of the chip-specific meta-data files or from \code{GOENTREZID2GO}.
}
\usage{
getOntology(inlist, ontology=c("MF", "BP", "CC"))
}
\arguments{
  \item{inlist}{ A list of GO associations}
  \item{ontology}{The name of the ontology you want returned. }
}
\details{
  The input list should be a list of lists, each element of \code{inlist}
  is itself a list containing the information that maps from a specified 
  ID (usually LocusLink) to GO information. Each element of the inner list 
  is a list with elements \code{GOID}, \code{Ontology} and \code{Evidence}. 
}
\value{
   A list of the same length as the input list. Each element of this
   list will contain a vector of \code{GOID}s for those terms that match
   the requested ontology.
}
\author{R. Gentleman}

\seealso{\code{\link{getEvidence}}, \code{\link{dropECode}} }

\examples{
 bb=GOENTREZID2GO$"4121"
 getOntology(bb)
 sapply(bb, function(x) x$Ontology)
}

}
\keyword{ manip }
