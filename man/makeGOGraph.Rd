\name{makeGOGraph}
\alias{makeGOGraph}
\title{ Construct a GO Graph }
\description{
  The directed acyclic graph (DAG) based on finding the most specific
  terms for the supplied LocusLink IDs is constructed and returned. The
  constructuion is per GO ontology (there are three, MF, BP and CC) and
  once the most specific terms have been identified then all less
  specific terms are found (these are the parents of the terms) and then
  their parents and so on, until the root is encountered.  
}
\usage{
makeGOGraph(x, Ontology="MF", removeRoot = TRUE)
}
\arguments{
  \item{x}{A vector of LocusLink IDs. }
  \item{Ontology}{Which of the GO ontologies to use. }
  \item{removeRoot}{A logical value indicating whether the GO root nood
  should be removed or not.}
}
\details{
  For each supplied LocusLink identifier all the GO annotations (in the
  specified ontology) are found.
  
  The mappings are different for the different ontologies. Typically a
  GO indentifier is used only in one specific ontology.

  The resulting structure is stored in a graph using the \code{graph}
  package, again from Bioconductor.
  
}
\value{
  An object that inherits from the \code{graph} class. The particular
  implementation is not specified.
}
\references{The Gene Ontology Consortium }
\author{R. Gentleman}

\seealso{\code{\link{oneGOGraph}}}

\examples{
 gN <- c("71733","224019", "237928", "11308")
 gg1 <- makeGOGraph(gN, "BP")
 
}
\keyword{manip}
