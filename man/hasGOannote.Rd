\name{hasGOannote}
\alias{hasGOannote}
\title{Check for GO annotation  }
\description{
  Given a GO term, or a vector of GO terms and an ontology this function
  determines which of the terms have GO annotation in the specified
  ontology. 

}
\usage{
hasGOannote(x, which="MF")
}
\arguments{
  \item{x}{A character vector, an instance of the \code{GOTerms} class
  or a list of \code{GOTerms}. }
  \item{which}{One of "MF", "BP" or "CC"}
}
\details{
  The available GO annotation is searched and a determination of
  whether a specific GO identifier has a value in the specified
  ontology is made.
}
\value{
  A logical vector of the same length as \code{x}.
}
\author{R. Gentleman}
\seealso{\code{\link{get}}}

\examples{
 t1 <- "GO:0003680"
 hasGOannote(t1)
 hasGOannote(t1, "BP")
}
\keyword{manip}
