\name{termGraphs}
\alias{termGraphs}
\alias{inducedTermGraph}
\alias{plotGOTermGraph}

\title{Extraction and Plotting of GO Terms from a GOHyperGResult Object}
\description{
  These functions extract and plot \code{graph} instances representing the
  relationships among GO terms tested using \code{hyperGTest}.
}
\usage{
termGraphs(r, id = NULL, pvalue = NULL, use.terms = TRUE)
inducedTermGraph(r, id, children = TRUE, parents = TRUE)
plotGOTermGraph(g, r = NULL, add.counts = TRUE, max.nchar = 20,
                node.colors=c(sig="lightgray", not="white"),
                node.shape="plaintext", ...)
}

\arguments{

  \item{r}{A \code{GOHyperGResult} object as returned by
    \code{hyperGTest} when given a \code{GOHyperGParams} object as
    input. }
  
  \item{id}{A character vector of category IDs that specifies which
    terms should be included in the graph.}

  \item{pvalue}{Numeric p-value cutoff to use for selecting category
    terms to include.  Will be ignored if \code{id} is present.}

  \item{use.terms}{Logical value indicating whether a \code{"term"} node
    attribute should be added to the returned graph providing the more
    descriptive, but possibly much longer, GO Terms.}

  \item{children}{A logical value indicating whether to include direct
    child terms of the terms specified by \code{id}.}

  \item{parents}{A logical value indicating whether to include direct
    parent terms of the terms specified by \code{id}.}

  \item{g}{A \code{graph} object as returned by \code{inducedTermGraph}
    or \code{termGraphs}.}

  \item{add.counts}{A logical value indicating whether category size
    counts should be added to the node labels when plotting.}

  \item{max.nchar}{The maximum character length for node labels in the plot.}

  \item{node.colors}{A named character vector of length two with
    compoents \code{sig} and \code{not}, giving color names for the
    significant and non-significant nodes, respectively.}

  \item{node.shape}{This argument controls the shape of the plotted
    nodes and must take on a value allowed by Rgraphviz.}

  \item{...}{For \code{plotGOTermGraph}, extra arguments are passed to
    the \code{plot} function.}

}

\details{
  \describe{
    
    \item{termGraphs}{returns a list of \code{graph} objects each
      representing one of the connected components of the subgraph of
      the GO ontology induced by selecting the specified GO IDs (if
      \code{id} is present) or by selecting the GO IDs that have a
      p-value less that \code{pvalue}.  If \code{use.terms} is
      \code{TRUE} the GO IDs will be translated into GO Term names and
      attached to the nodes as node attributes (see \code{nodeData}).
      Edges in the graphs go from child (more specific) to parent (less
      specific).}

    \item{inducedTermGraph}{returns a \code{graph} object representing
      the GO graph induced by the terms specified by \code{id}.  The
      \code{children} and \code{parent} arguments control whether direct
    children and/or direct parents of the terms specified by \code{id}
    are added to the graph (at least one of the two must be
    \code{TRUE}).}

    \item{plotGOTermGraph}{Create a plot using Rgraphviz of a
      \code{graph} object as returned by either \code{termGraphs} or
      \code{inducedTermGraph}.  If a \code{GOHyperGResult} object is
      provided, then the nodes will be colored according to significance
      (based on the result object's \code{pvalueCutoff}) and counts will
      be added to show the size of the categories.
    }

  }
}

\author{Seth Falcon}

\keyword{ manip }

