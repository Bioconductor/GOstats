\name{oneGOGraph}
\alias{oneGOGraph}
\title{Build a GO graph for one identifier  }
\description{
  Given a single GO identifier and a set of mappings
 to the less specific sets of nodes this function will construct 
 the graph that includes that node and all children down to the root node
 for the ontology.
}
\usage{
oneGOGraph(x, dataenv)
}
\arguments{
  \item{x}{A length one character vector with the name of the term. }
  \item{dataenv}{ An environment for finding the parents of that term. }
}
\details{
  For any gene we define the induced GO graph to be that graph, based on
  the DAG structure (child - parent) of the GO ontology of terms.
}
\value{
  The induced GO graph (or NULL) for the given GO identifier.
}
\author{R. Gentleman}

\seealso{\code{\link{makeGOGraph}}}

\examples{

 g1 <- oneGOGraph("GO:0003680", GOMFPARENTS)
 g2 <- oneGOGraph("GO:0003701", GOMFPARENTS)
 g3 <- combGOGraph(g1, g2)

if( require(Rgraphviz) && interactive() )
  plot(g3)
}
\keyword{manip}
