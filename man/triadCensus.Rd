\name{triadCensus}
\alias{triadCensus}
\alias{triad}
\alias{enumPairs}
\alias{isTriad}
\alias{reduce2Degreek}
\title{ Triad Functions }
\description{
  These functions provide some tools for finding triads in an undirected
  graph. A triad is a clique of size 3. The function \code{triadCensus}
  returns a list of all triads.
}
\usage{
triadCensus(graph)
isTriad(x, y, z, elz, ely)
reduce2Degreek(graph, k)
enumPairs(iVec)

}
\arguments{
  \item{graph}{An instance of the \code{graph} class. }
  \item{k}{An integer indicating the minimum degree wanted.}
  \item{x}{A node}
  \item{y}{A node}
  \item{z}{A node}
  \item{elz}{The edgelist for \code{z}}
  \item{ely}{The edgelist for \code{y}}
  \item{iVec}{A vector of unique values}
}
\details{
  \code{enumPairs} takes a vector as input and returns a list of length
  \code{choose(length(iVec),2)/2} containing all unordered pairs of
  elements.

  \code{isTriad} takes three nodes as arguments. It is already known
  that \code{x} has edges to both \code{y} and \code{z} and we want to
  determine whether these are reciprocated. This is determined by
  examining \code{elz} for both \code{x} and \code{y} and then examining
  \code{ely} for both \code{x} and \code{z}. 

  \code{reduce2Degreek} is a function that takes an undirected graph as
  input and removes all nodes of degree less than \code{k}. This process
  is iterated until there are no nodes left (an error is thrown) or all
  nodes remaining have degree at least \code{k}. The resultant subgraph
  is returned. It is used here because to be in a triad all nodes must
  have degree 2 or more.

  \code{triadCensus} makes use of the helper functions described above
  and finds all triads in the graph.
  
}
\value{
  A list where each element is a triple indicating the members of the
  triad. Order is not important and all triads are reported in
  alphabetic order.
}

\author{R. Gentleman}
\note{See the graph package, RBGL and Rgraphviz for more details and
  alternatives. }


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
}
\keyword{manip}
