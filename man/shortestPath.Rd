\name{shortestPath}
\alias{shortestPath}
\title{ Shortest Path Analysis }
\description{
 The shortest path analysis was proposed by Zhou et. al. The basic
 computation is to find the shortest path in a supplied graph between
 two Entrez Gene IDs. Zhou et al claim that other genes annotated along
 that path are likely to have the same GO annotation as the two end
 points.
}
\usage{
shortestPath(g, GOnode, mapfun=NULL, chip=NULL)
}
\arguments{
  \item{g}{An instance of the \code{graph} class. }
  \item{GOnode}{A length one character vector specifying the GO node of
    interest. }

    \item{mapfun}{A function taking a character vector of GO IDs as its
    only argument and returning a list of character vectors of Enterz
    Gene IDs annotated at each corresponding GO ID.  The function should
    behave similarly to \code{mget(x, go2egmap, ifnotfound=NA)}, that
    is, \code{NA} should be returned if a specified GO ID has no Entrez
    ID mappings.  See details for the interaction of \code{mapfun} and
    \code{chip}.}

  \item{chip}{The name of a DB-based annotation data package (the name
    will end in ".db").  This package will be used to generate an Entrez
    ID to GO ID mapping instead of \code{mapfun}.}
}
\details{
  The algorithm implemented here is quite simple. All Entrez Gene 
  identifiers that are annotated at the GO node of interest are
  obtained. Those that are found as nodes in the graph are retained and
  used for the computation. For every pair of nodes at the GO term the
  shortest path between them is computed using \code{sp.between} from
  the RBGL package.

  There is a presumption that the graph is \code{undirected}. This
  restriction could probably be lifted if there was some reason for it -
  a patch would be gratefully accepted.

  The mapping of GO node to Entrez ID is achieved in one of three ways:

  \enumerate{
    \item If \code{mapfun} is provided, it will be used to perform the
      needed lookups.  In this case, \code{chip} will be ignored.

    \item If \code{chip} is provided and \code{mapfun=NULL}, then the
      needed lookups will be done based on the GO to Entrez mappings
      encapsulated in the specified annotation data package.  This is
      the recommended usage.

    \item If \code{mapfun} and \code{chip} are \code{NULL} or missing,
      then the function will attempt to load the GO package (the
      environment-based package, distinct from GO.db).  This package
      contains a legacy environment mapping GO IDs to Entrez IDs.  If
      the GO package is not available, an error will be raised.
      Omitting both \code{mapfun} and \code{chip} is not recommended as
      it is not compatible with the DB-based annotation data packages.
  }
  
}
\value{
  The return values is a list with the following components:
  \item{shortestpaths }{A list of the ouput from \code{sp.between}. The
    names are the names of the nodes used as the two endpoints}
  \item{nodesUsed }{A vector of the Entrez Gene IDs that were both found
    at the GO term of interest and were nodes in the supplied graph,
    \code{g}. These were used to compute the shortest paths.}
  \item{nodesNotUsed}{A vector of Entrez Gene IDs that were annotated at
    the GO term, but were not found in the graph \code{g}.}
}
\references{Transitive functional annotation by shortest-path analysis
  of gene expression data, by X. Zhou and M-C J. Kao and W. H. Wong,
  PNAS, 2002}

\author{R. Gentleman }

\seealso{\code{\link[RBGL]{sp.between}}}

\examples{
library("hgu95av2.db")
library("RBGL")

set.seed(321)
uniqun <- function(x) unique(unlist(x))

goid <- "GO:0005778"
egIds <- uniqun(mget(uniqun(hgu95av2GO2PROBE[[goid]]),
                            hgu95av2ENTREZID))

v1 <- randomGraph(egIds, 1:10, .3, weights=FALSE)
## Since v1 is random, it might be disconnected and we need a
## connected graph to guarantee the existence of a path.
c1 <- connComp(v1)
largestComp <- c1[[which.max(sapply(c1, length))]]
v2 <- subGraph(largestComp, v1)

a1 <- shortestPath(v2, goid, chip="hgu95av2.db")

}
\keyword{ manip }
