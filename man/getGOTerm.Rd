\name{getGOTerm}
\alias{getGOTerm}
\alias{getGOParents}
\alias{getGOChildren}
\alias{getGOCategory}

\title{Functions to Access GO data.  }
\description{
  These functions provide access to data in the GO package. The data are
  assembled from publically available data from the Gene Ontology
  Consortium (GO), \url{www.go.org}. Given a list of GO identifiers they
  access the children (more specific terms), the parents (less specific
  terms) and the terms themselves.
}
\usage{
getGOTerm(x)
getGOParents(x)
getGOChildren(x)
getGOCategory(x)
}
\arguments{
  \item{x}{A character vector of valid GO identifiers. }
}
\details{
  GO consists of three (soon to be more) specific hierarchies: Molecular
  Function (MF), Biological Process (BP) and Cellular Component
  (CC). For more details consult the GO website. For each GO identifier
  each of these three hierarchies is searched and depending on the
  function called the appropriate values are obtained and returned.

  It is possible for a GO identifier to have no children or for it to
  have no parents. However, it must have a term associated with it.
}
\value{
  A list of the same length as \code{x}.
  The list contains one entry for each element of \code{x}. That entry
  is itself a list. With one component named \code{Ontology} which
  has as its value one of MF, BP or CC. The second component has the
  same name as the suffix of the call, i.e. Children, Parents, or Term.
  If there was no match in any of the ontologies then a length zero list
  is returned for that element of \code{x}.

  For \code{getGOCategory} a vector of categories (the names of which
  are the original GO term names). Elements of this list that are
  \code{NA} indicate term names for which there is no category (and
  hence they are not really term names).
}
\references{The Gene Ontology Consortium}
\author{R. Gentleman}


\examples{
 library(GO)

 sG <- c("GO:0005515", "GO:0000123", "GO:0000124", "GO:0000125",
             "GO:0000126", "GO:0020033", "GO:0006830", "GO:0009874",
             "GO:0015916", "GO:0015339")


 gT <- getGOTerm(sG)

 gP <- getGOParents(sG)

 gC <- getGOChildren(sG)

 gcat <- getGOCategory(sG)

}
\keyword{manip}
