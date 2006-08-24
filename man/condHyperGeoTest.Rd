\name{condHyperGeoTest}
\alias{condHyperGeoTest}
\alias{condHyperGeoTest,GeneGoCondHyperGeoTestParams-method}

\title{Conditional Hypergeometric Tests for GO}

\description{
  Given a set of unique Entrez Gene Identifiers, a microarray
  annotation data package name, and the GO ontology of interest, this
  function will compute conditional Hypergeomtric p-values for
  overrepresentation of each GO term in the specified ontology among
  the GO annotations for the interesting genes (as indicated by the
  Entrez Gene ids).
}
\usage{
condHyperGeoTest(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{A \code{GeneGoHyperGeoTestParams} instance }
  \item{\dots}{ unused }
  \item{pval.cutoff}{Hypergeometric p-values smaller than
      \code{pval.cutoff} will be deemed significant.  }

}
\details{
  FIXME
}
\value{
  A \code{GeneGoCondHyperGeoTestResult} instance.
}
\references{ FIXME}
\author{Seth Falcon}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\keyword{models}

