\name{makeGOGraph}
\alias{makeGOGraph}
\title{ Construct a GO Graph }
\description{
  The directed acyclic graph (DAG) based on finding the most specific
  terms for the supplied identifiers is constructed and returned. The
  constructuion is per GO ontology (there are three, MF, BP and CC) and
  once the most specific terms have been identified then all less
  specific terms are found (these are the parents of the terms) and then
  their parents and so on, until the root is encountered.  
}
\usage{
makeGOGraph(x, what="MF", lib="hgu95av2", removeRoot = TRUE)
}
\arguments{
  \item{x}{A vector of identifiers. }
  \item{what}{Which of the GO ontologies to use. }
  \item{lib}{The name of a meta-data package to use for mapping the
    identifiers to GO identifiers. }
  \item{removeRoot}{A logical value indicating whether the GO root nood
  should be removed or not.}
}
\details{
  The mapping of manufacturers identifiers (e.g. Affymetrix) to GO
  identifiers is done on the basis of an initial mapping to LocusLink
  identifiers. For many data sources these are available from the
  Bioconductor Project (\url{www.bioconductor.org}).

  Once that mapping has occurred we obtain from the meta-data package
  the mapping to the most specific terms (again these have been
  precomputed) and using those together with the GO package (again from
  Bioconductor) mapping to all parents down to the root node is
  performed.

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
 gN <- c("38940_at","2073_s_at", "35580_at",  "34701_at")
 gg1 <- makeGOGraph(gN, "BP", "hgu95av2")
 
}
\keyword{manip}
