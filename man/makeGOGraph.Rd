\name{makeGOGraph}
\alias{makeGOGraph}
\title{ Construct a GO Graph }
\description{
  The directed acyclic graph (DAG) based on finding the most specific
  terms for the supplied Entrez Gene IDs is constructed and returned. The
  constructuion is per GO ontology (there are three, MF, BP and CC) and
  once the most specific terms have been identified then all less
  specific terms are found (these are the parents of the terms) and then
  their parents and so on, until the root is encountered.  
}
\usage{
makeGOGraph(x, Ontology, eg2go.map, removeRoot = TRUE)
}
\arguments{
  \item{x}{A vector of Entrez Gene IDs. }
  \item{Ontology}{Which GO ontology to use (CC, BP, or MF). }
  \item{eg2go.map}{An object with an environment-like interface
    with an \code{mget} method that provides a mapping from Entrez Gene
    IDs to GO IDs.}
  \item{removeRoot}{A logical value indicating whether the GO root node
  should be removed or not.}
}
\details{
  For each supplied Entrez Gene identifier all the GO annotations (in the
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
 library("hgu95av2.db")
 set.seed(321)
 uniqun <- function(x) unique(unlist(x))
 gN <- uniqun(sample(hgu95av2ENTREZID, 4))
 probes <- uniqun(mget(gN, revmap(hgu95av2ENTREZID)))
 goids <- mget(probes, hgu95av2GO)
 names(goids) <- gN
 eg2go <- l2e(goids)
 gg1 <- makeGOGraph(gN, "BP", eg2go)
 
}
\keyword{manip}
