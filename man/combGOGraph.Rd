\name{combGOGraph}
\alias{combGOGraph}
\title{Combine GO Graphs  }
\description{
  Given two GO graphs this function combines them into a single graph
  and returns that.
}
\usage{
combGOGraph(g1, g2)
}
\arguments{
  \item{g1}{An instance of a \code{graph}.}
  \item{g2}{An instance of a \code{graph}. }
}
\details{
  A new \code{graph} object is created by linking the two supplied
  graphs together. Since they are both GO graphs they will share, at
  least, the GO root node.

  This could probably be replaced by the \code{union} function in the
  \code{graph} package.
}
\value{
  A graph containing the union of the nodes and edges in the two
  graphs.
}
\author{R. Gentleman}

\seealso{\code{\link{union}}}

\examples{
 g1 <- oneGOGraph("GO:0003680", GOMFPARENTS)
 g2 <- oneGOGraph("GO:0003701", GOMFPARENTS)
 g3 <- combGOGraph(g1, g2)
}
\keyword{manip}
