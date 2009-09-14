\name{hyperGTest}
\alias{hyperGTest}
\alias{hyperGTest,GOHyperGParams-method}

\title{Hypergeometric Tests for GO term association}

\description{ Given a \code{GOHyperGParams} instance
  containing a set of unique Entrez Gene Identifiers, a microarray
  annotation data package name, and the GO ontology of interest, this
  function will compute Hypergeomtric p-values for over or
  under-representation of each GO term in the specified ontology among
  the GO annotations for the interesting genes.  The computations can
  be done conditionally based on the structure of the GO graph.}


\arguments{
  \item{p}{A \code{GOHyperGParams} instance }
}

\details{

  When \code{conditional(p) == TRUE}, the \code{hyperGTest} function
  uses the structure of the GO graph to estimate for each term whether
  or not there is evidence beyond that which is provided by the term's
  children to call the term in question statistically overrepresented.

  The algorithm conditions on all child terms that are themselves
  significant at the specified p-value cutoff.  Given a subgraph of
  one of the three GO ontologies, the terms with no child categories
  are tested first.  Next the nodes whose children have already been
  tested are tested.  If any of a given node's children tested
  significant, the appropriate conditioning is performed.
}
\value{
  A \code{GOHyperGResult} instance.
}

\references{FIXME}

\author{Seth Falcon}
\seealso{
  \code{\link{GOHyperGResult-class}},
  %\code{\link[Category]{geneCategoryHyperGeoTest}},
  \code{\link[Category]{geneGoHyperGeoTest}},
  \code{\link[Category]{geneKeggHyperGeoTest}}
}
\keyword{models}

