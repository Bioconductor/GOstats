\name{simLL}
\alias{simLL}
\alias{simUI}
\alias{simLP}
\title{Functions to compute similarities between GO graphs and also
between LocusLink IDs based on their induced GO graphs.}
\description{
  Both \code{simUI} and \code{simLP} compute a similarity measure
  between two GO graphs. For \code{simLL}, first the induced GO graph
  for each of its arguments is found and then these are passed to one
  of \code{simUI} or \code{simLP}.
}
\usage{
simLL(ll1, ll2, Ontology = "MF", measure = "LP", dropCodes = NULL)
simUI(g1, g2)
simLP(g1, g2)
}

\arguments{
  \item{ll1}{A LocusLink ID as a character vector. }
  \item{ll2}{A LocusLink ID as a character vector.}
  \item{Ontology}{Which ontology to use ("MF", "BP", "CC"). }
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

  For code{\simLP} the length of the longest path in the intersection
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
}

\author{R. Gentleman }

\seealso{\code{\link{makeGOGraph}} }
\examples{
  aa = simLL("9184", "3547", "BP")
}
\keyword{manip}

