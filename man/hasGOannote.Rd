\name{hasGOannote}
\alias{hasGOannote}
\title{Check for GO annotation  }
\description{
  Given a GO term and an environment this function determines whether
  there is an entry (or symbol) in the environment with the same name as
  the term.
}
\usage{
hasGOannote(x, which="MF")
}
\arguments{
  \item{x}{A length one character vector. }
  \item{which}{Either an environment or the name of an environment. By
    default this is the MF parents environment. }
}
\details{
  The element \code{x} is searched for in the environment. If it is
  found then \code{TRUE} is returned otherwise \code{FALSE} is returned.
}
\value{
  As described above.
}
\author{R. Gentleman}
\seealso{\code{\link{get}}}

\examples{
 t1 <- "GO:0003680"
 hasGOannote(t1)
 hasGOannote(t1, "BP")
}
\keyword{manip}
