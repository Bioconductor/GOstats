\name{GOHyperG}
\alias{GOHyperG}
\title{ Hypergeometric Tests for GO  }
\description{
 Given a set of unique LocusLink Identifiers, a microarray chip and the
 GO category of interest this function will compute all Hypergeomtric
 p-values for overrepresentation of the interesting genes (as indicated
 by the unique LocusLink Identifiers) at the nodes in the induced GO
 graph.
}
\usage{
GOHyperG(x, lib="hgu95av2", what="MF")
}
\arguments{
  \item{x}{A vector of unique LocusLink identifiers. }
  \item{lib}{The name of the annotation library for the chip that was used. }
  \item{what}{One of "MF", "BP", or "CC" indicating which of the GO
    categories the computations should be made for.}
}
\details{
  Typical usage will be to have a microarray experiment from which a set
  of interesting genes/probes has been obtained. To determine whether
  there is an overrepresentation of these genes at particular GO terms
  a simple hypergeometric calculation has often been made. Two
  substantial issues arise. First and most importantly it is not clear
  how to do any form of p-value correction in this case. The tests are
  not independent and the underlying structure of the GO graph present
  certain problems that still need to be addressed. The second substantial
  issue is that arises is that the mappings are based on LocusLink
  identifiers and hence all computations should also be based on
  unique LocusLink identifiers. In \code{GOHyperG} every attempt to
  appropriately correct for non-uniqueness of mappings has been made.
  
  The user provides a vector of unique LocusLink identifiers and these
  are used, together with the name of the chip to create the necessary
  counts. It is important that the correct chip be identified as that
  determines the overall counts and all inference will be incorrect if
  that is not correct.

  The test performed is a Hypergeometric test, using \code{phyper},
  where at each GO node we determine how many LLIDs from the chip were
  annotated there, how many of the supplied LLIDs were annotated there
  and compute a $p$-value. This is the equivalent of using Fisher's
  exact test.
  
}
\value{
  The returned value is a list with components:
  \item{pvalues }{The ordered p-values.}
  \item{goCounts}{The vector of counts of LLIDs from the chip at each node.}
  \item{intCounts}{The vector of counts of the supplied LLIDs annotated
    at each node.}
  \item{numLL}{The number of unique LLIDs on the chip that are mapped to
    some term in the specified GO category.}
  \item{numInt}{The number of unique LLIDs from those supplied that are
    mapped to some term in the specified GO category.}
}
\author{R. Gentleman}

\seealso{\code{\link{phyper}}}

\examples{

library(hgu95av2)
library(GO)
w1<-as.list(hgu95av2LOCUSID)
w2<-unique(unlist(w1))
set.seed(123)
#pick a hundred interesting genes
 myLL <- sample(w2, 100)
 xx<-GOHyperG(myLL)
xx$numLL
xx$numInt
sum(xx$pvalues < 0.01)

}
\keyword{ htest }
