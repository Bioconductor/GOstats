\name{compGdist}
\alias{compGdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to compute the distance between pairs of nodes in a graph. }
\description{
  Given a graph, \code{g}, and a set of nodes in the graph,
  \code{whNodes}, Dijkstra's shortest path algorithm is used to compute
  the distance between all pairs of nodes in \code{whNodes}.
}
\usage{
compGdist(g, whNodes, verbose = FALSE)
}
\arguments{
  \item{g}{ An instance of the \code{graph} class. }
  \item{whNodes}{A vector of lables of the nodes in \code{g} for which
  distances are to be computed. }
  \item{verbose}{If \code{TRUE} then output reporting the progress will
  be reported. }
}
\details{
  This function can be quite slow, computation of the pairwise
  distances is not especially fast and if \code{whNodes} is long then
  there are many of them to compute.
}
\value{
  A matrix containing the pairwise distances. It might be worth making
  this an instance of the \code{dist} class at some point.
}

\author{R. Gentleman }

\seealso{ \code{\link{compCorrGraph}} }
\examples{

 ##none yet
}
\keyword{manip }

