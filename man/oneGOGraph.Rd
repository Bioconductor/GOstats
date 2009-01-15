\name{oneGOGraph}
\alias{oneGOGraph}
\alias{GOGraph}
\alias{getGoGraph}
\alias{getGoGraph,GOHyperGParams,character-method}

\title{Construct the GO graph given a set of leaves.  }

\description{
  Given one or more GO identifiers (which indicate the leaves in the
  graph) and a set of mappings to the less specific sets of nodes this
  function will construct the graph that includes that node and all
  children down to the root node for the ontology.
}
\usage{
oneGOGraph(x, dataenv)
GOGraph(x, dataenv)
getGoGraph(p, goIds)

}
\arguments{
  \item{x}{A character vector of GO identifiers. }
  \item{dataenv}{ An environment for finding the parents of that
      term. }

  \item{p}{A \code{GOHyperGParams} instance.}
  \item{goIds}{A character vector providing the GO term IDs.}

}
\details{
  For any set of GO identifiers (from a common ontology) we define the
  induced GO graph to be  that graph, based on 
  the DAG structure (child - parent) of the GO ontology of terms, which
  takes the most specific set of GO terms that apply (for that
  ontology) and then joins these to all less specific terms. These
  functions help construct such graphs.
}
\value{
  The induced GO graph (or NULL) for the given GO identifier.
}
\author{R. Gentleman}

\seealso{\code{\link{makeGOGraph}}}

\examples{

 library("GO.db")
 g1 <- oneGOGraph("GO:0003680", GOMFPARENTS)
 g2 <- oneGOGraph("GO:0003701", GOMFPARENTS)
 g3 <- join(g1, g2)

 g4 <- GOGraph(c("GO:0003680", "GO:0003701"), GOMFPARENTS)
if( require("Rgraphviz") && interactive() )
  plot(g3)
}
\keyword{manip}
