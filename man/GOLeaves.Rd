\name{GOLeaves}
\alias{GOLeaves}
\title{Identify the leaves in a GO Graph  }
\description{
  Given a GO graph this function returns the node labels for all leaves
  in the graph. A leaf is defined to be a node with only out-edges and
  no in-edges.
}
\usage{
GOLeaves(inG)
}
\arguments{
  \item{inG}{An instance of a GO graph. }
}
\details{
  All nodes in \code{inG} are inspected for in-edges and those with none
  are returned.

  This should probably be replaced by a function in the graph package
  that identifies leaves, there is nothing special about GO here.
  
}
\value{
  A vector of the node labels for those nodes with no in-edges.
}

\author{R. Gentleman}

\seealso{\code{\link{makeGOGraph}}}

\examples{
 g1 <- oneGOGraph("GO:0003680", GOMFPARENTS)
 g2 <- oneGOGraph("GO:0003701", GOMFPARENTS)
 g3 <- combGOGraph(g1, g2)
 GOLeaves(g3)
}
\keyword{manip}
