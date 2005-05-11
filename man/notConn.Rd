\name{notConn}
\alias{notConn}
\title{Find genes that are not connected to the others. }
\description{
  A function that takes as input a distance matrix and finds those
  entries that are not connected to any others (ie. those with distance
  \code{Inf}. 
}
\usage{
notConn(dists)
}
\arguments{
  \item{dists}{A distance matrix.}
}
\details{
  It is a very naive implementation. It presumes that not connected
  entries are not connected to any other entries, and this might not be
  true. Using the \code{connComp} function from the \code{graph}
  package or the \code{RBGL} package might be a better approach.
}
\value{
 A vector of the names of the items that are not connected.
}
\author{R. Gentleman }
\seealso{\code{\link[graph]{connComp}}}
\examples{
 data(Ndists)
 notConn(Ndists)
}
\keyword{manip}
