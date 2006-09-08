\name{goHyperGeoTest}
\alias{goHyperGeoTest}
\alias{goHyperGeoTest,GeneGoCondHyperGeoTestParams-method}

\title{Conditional Hypergeometric Tests for GO}

\description{ Given a \code{GeneGoCondHyperGeoTestParams} instance
  containing a set of unique Entrez Gene Identifiers, a microarray
  annotation data package name, and the GO ontology of interest, this
  function will compute conditional Hypergeomtric p-values for
  overrepresentation of each GO term in the specified ontology among
  the GO annotations for the interesting genes.} 

\usage{ goHyperGeoTest(p) }

\arguments{
  \item{p}{A \code{GeneGoHyperGeoTestParams} instance }
}

\details{

  The \code{goHyperGeoTest} function uses the structure of the GO
  graph to estimate for each term whether or not there is evidence
  beyond that which is provided by the term's children to call the
  term in question statistically overrepresented.

  The algorithm conditions on all child terms that are themselves
  significant at a given p-value cutoff.  Given a subgraph of one of
  the three GO ontologies, the terms with no child categories are
  tested first.  Next the nodes whose children have already been
  tested are tested.  If any of a given node's children tested
  significant, the appropriate conditioning is performed.
}
\value{
  A \code{GeneGoHyperGeoTestResult} instance.
}

\references{FIXME}

\author{Seth Falcon}
\seealso{
  \code{\link{GeneGoHyperGeoTestResult-class}}
  \code{\link[Category]{geneCategoryHyperGeoTest}}
  \code{\link[Category]{geneGoHyperGeoTest}}
  \code{\link[Category]{geneKeggHyperGeoTest}}
}
\keyword{models}

