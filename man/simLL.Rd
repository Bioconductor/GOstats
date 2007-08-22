\name{simLL}
\alias{simLL}
\alias{simUI}
\alias{simLP}
\title{Functions to compute similarities between GO graphs and also
between Entrez Gene IDs based on their induced GO graphs.}
\description{
  Both \code{simUI} and \code{simLP} compute a similarity measure
  between two GO graphs. For \code{simLL}, first the induced GO graph
  for each of its arguments is found and then these are passed to one
  of \code{simUI} or \code{simLP}.
}
\usage{
simLL(ll1, ll2, Ontology = "MF", measure = "LP", dropCodes = NULL,
      mapfun = NULL, chip = NULL)
simUI(g1, g2)
simLP(g1, g2)
}

\arguments{
  \item{ll1}{A Entrez Gene ID as a character vector. }
  \item{ll2}{A Entrez Gene ID as a character vector.}
  \item{Ontology}{Which ontology to use ("MF", "BP", "CC"). }
  \item{measure}{ Which measure to use ("LP", "UI"). }
  \item{dropCodes}{ A set of evidence codes to be ignored in
    constructing the induced GO graphs. }

  \item{mapfun}{A function taking a character vector of Entrez Gene IDs
    as its only argument and returning a list of "GO lists" matching the
    structure of the lists in the GO maps of annotation data packages.
    The function should behave similarly to \code{mget(x, eg2gomap,
      ifnotfound=NA)}, that is, \code{NA} should be returned if a
      specified
    Entrez ID has no GO mapping.  See details for the interaction of
    \code{mapfun} and \code{chip}.}

  \item{chip}{The name of a DB-based annotation data package (the name
    will end in ".db").  This package will be used to generate an Entrez
    ID to GO ID mapping instead of \code{mapfun}.}

  \item{g1}{An instance of the \code{graph} class.}
  \item{g2}{An instance of the \code{graph} class.}
}
\details{
  For each of \code{ll1} and \code{ll2} the set of most specific GO
  terms within the ontology specified (\code{Ontology}) that are not
  based on any excluded evidence code (\code{dropCodes}) are found.  The
  mapping is achieved in one of three ways:

  \enumerate{
    \item If \code{mapfun} is provided, it will be used to perform the
      needed lookups.  In this case, \code{chip} will be ignored.

    \item If \code{chip} is provided and \code{mapfun=NULL}, then the
      needed lookups will be done based on the Entrez to GO mappings
      encapsulated in the specified annotation data package.  This is
      the recommended usage.

    \item If \code{mapfun} and \code{chip} are \code{NULL} or missing,
      then the function will attempt to load the GO package (the
      environment-based package, distinct from GO.db).  This package
      contains a legacy environment mapping Entrez IDs to GO IDs.  If
      the GO package is not available, an error will be raised.
      Omitting both \code{mapfun} and \code{chip} is not recommended as
      it is not compatible with the DB-based annotation data packages.
  }

  Next, the induced GO graphs are computed.

  Finally these graphs are passed to one of \code{simUI}, (union
  intersection), or \code{simLP} (longest path). For \code{simUI} the
  distance is the size of the intersection of the node sets divided by
  the size of the union of the node sets. Large values indicate more
  similarity. These similarities are between 0 and 1.

  For \code{simLP} the length of the longest path in the intersection
  graph of the two supplied graph. Again, large values indicate more
  similarity. Similarities are between 0 and the maximum leaf depth of
  the graph for the specified ontology.
}
\value{
  A list with:
  \item{sim }{The numeric similarity measure.}
  \item{measure}{Which measure was used.}
  \item{g1}{The graph induced by \code{ll1}.}
  \item{g2}{The graph induced by \code{ll2}.}

  If one of the supplied Gene IDs does not have any GO terms associated
with it, in the selected ontology and with the selected evidence codes
then \code{NA} is returned.
}

\author{R. Gentleman }

\seealso{\code{\link{makeGOGraph}} }
\examples{
  library("hgu95av2.db")
  eg1 = c("9184", "3547")

  bb = simLL(eg1[1], eg1[2], "BP", chip="hgu95av2.db")
}
\keyword{manip}

