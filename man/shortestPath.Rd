\name{shortestPath}
\alias{shortestPath}
\title{ Shortest Path Analysis }
\description{
 The shortest path analysis was proposed by Zhou et. al. The basic
 computation is to find the shortest path in a supplied graph between
 two LocusLink IDs. Zhou et al claim that other genes annotated along
 that path are likely to have the same GO annotation as the two end
 points.
}
\usage{
shortestPath(g, GOnode)
}
\arguments{
  \item{g}{An instance of the \code{graph} class. }
  \item{GOnode}{A length one character vector specifying the GO node of
    interest. }
}
\details{
  The algorithm implemented here is quite simple. All LocusLink
  identifiers that are annotated at the GO node of interest are
  obtained. Those that are found as nodes in the graph are retained and
  used for the computation. For every pair of nodes at the GO term the
  shortest path between them is computed using \code{sp.between} from
  the RBGL package.

  There is a presumption that the graph is \code{undirected}. This
  restriction could probably be lifted if there was some reason for it -
  a patch would be gratefully accepted.
  
}
\value{
  The return values is a list with the following components:
  \item{shortestpaths }{A list of the ouput from \code{sp.between}. The
    names are the names of the nodes used as the two endpoints}
  \item{nodesUsed }{A vector of the LocusLink IDs that were both found
    at the GO term of interest and were nodes in the supplied graph,
    \code{g}. These were used to compute the shortest paths.}
  \item{nodesNotUsed}{A vector of LocusLink IDs that were annotated at
    the GO term, but were not found in the graph \code{g}.}
}
\references{Transitive functional annotation by shortest-path analysis
  of gene expression data, by X. Zhou and M-C J. Kao and W. H. Wong,
  PNAS, 2002}

\author{R. Gentleman }

\seealso{\code{\link{sp.between}}}

\examples{
library(GO)
library(RBGL)

tst <- unique(unlist(mget(c("GO:0005778", "GO:0005779",
                            "GO:0030060"), GOGO2LL)))
set.seed(123)
v1 <- randomGraph(tst, 1:10, .3)


a1 <- shortestPath(v1, "GO:0005779")


}
\keyword{ manip }
