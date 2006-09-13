\name{compCorrGraph}
\alias{compCorrGraph}
\title{A function to compute a correlation based graph from Gene
Expression Data }
\description{
  Given a set of gene expression data (an instance of the
  \code{exprSet} class) this function computes a graph based on
  correlations between the probes.
}
\usage{
compCorrGraph(eSet, k = 1, tau = 0.6)
}
\arguments{
  \item{eSet}{An instance of the \code{exprSet} class. }
  \item{k}{The power to raise the correlations to. }
  \item{tau}{The lower cutoff for absolute correlations. }
}
\details{
  Zhou et al. describe a method of computing a graph between probes
  (genes) based on estimated correlations between probes. This function
  implements some of their methods.

  Pearson correlations between probes are computed and then these are
  raised to the power \code{k}. Any of the resulting estimates that are
  less than \code{tau} in absolute value are set to zero. 
}
\value{
  An instance of the \code{graph} class. With edges and edge weights
  determined by applying the algorithm described previously.
}
\references{Zhou et al., Transitive functional annotation by
shortest-path analysis of gene expression data.}
\author{R. Gentleman}

\seealso{\code{\link{compGdist}}}
\examples{

 ## Create an ExpressionSet to work with
 set.seed(123)
 exprMat <- matrix(runif(100 * 10), nrow=100)
 genData <- new("ExpressionSet", exprs=exprMat)

 corrG = compCorrGraph(genData)

}
\keyword{ manip }

