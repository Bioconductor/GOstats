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
simLL(ll1, ll2, Ontology, eg2go.map, measure = "LP", dropCodes = NULL)
simUI(g1, g2)
simLP(g1, g2)
}

\arguments{
  \item{ll1}{A Entrez Gene ID as a character vector. }
  \item{ll2}{A Entrez Gene ID as a character vector.}
  \item{Ontology}{Which ontology to use ("MF", "BP", "CC"). }
  \item{eg2go.map}{An object with an environment-like interface
    with an \code{mget} method that provides a mapping from Entrez Gene
    IDs to GO IDs.}
  \item{measure}{ Which measure to use ("LP", "UI"). }
  \item{dropCodes}{ A set of evidence codes to be ignored in
    constructing the induced GO graphs. }
  \item{g1}{An instance of the \code{graph} class.}
  \item{g2}{An instance of the \code{graph} class.}
}
\details{
  For each of \code{ll1} and \code{ll2} the set of most specific GO
  terms within the ontology specified (\code{Ontology}) that are not
  based on any excluded evidence code (\code{dropCodes}) are
  found. Then the induced GO graphs are computed.

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

  probes <- unique(unlist(mget(eg1, revmap(hgu95av2ENTREZID))))
  goids <- mget(probes, hgu95av2GO)
  names(goids) <- eg1
  eg2go <- l2e(goids)

  bb = simLL(eg1[1], eg1[2], "BP", eg2go)
}
\keyword{manip}

