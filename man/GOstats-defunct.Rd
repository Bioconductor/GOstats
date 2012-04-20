\name{GOstats-defunct}
\alias{GOstats-defunct}
\title{Defunct Functions in GOstats Package}
%------ PLEASE: put \alias{.} here for  EACH !
\alias{combGOGraph}
\alias{hyperGtable}
\alias{hyperG2Affy}
\alias{GOLeaves}
\alias{selectedGenes}
\alias{GOHyperG}
\alias{GOKEGGHyperG}
\alias{getGoGraph}

\description{
  The functions or variables listed here are no longer part of GOstats as
  they are not needed (any more).
}
\usage{
combGOGraph()
hyperGtable()
hyperG2Affy()
selectedGenes()
GOHyperG()
GOKEGGHyperG()
getGoGraph()
}
\details{
  \code{combGOGraph} was replaced by \code{join}.
  \code{hyperGtable} was replaced by \code{summary}.
  \code{hyperG2Affy} was replaced by \code{probeSetSummary}.
  \code{GOLeaves} was replaced by \code{graph::leaves}.
  \code{selectedGenes} was replaced by \code{geneIdsByCategory}.
  \code{GOHyperG} was replaced by \code{hyperGTest}.
  \code{GOKEGGHyperG} was replaced by \code{hyperGTest}.
  \code{getGoGraph} was replaced by \code{GOGraph}.
}

\keyword{manip}
