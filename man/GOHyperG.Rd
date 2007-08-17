\name{GOHyperG}
\alias{GOHyperG}

\title{(DEPRECATED) Hypergeometric Tests for GO  }
\description{
  Use \code{hyperGTest} instead.
  
 Given a set of unique Entrez Gene Identifiers, a microarray annotation
 data package name, and the GO category of interest, this function will
 compute Hypergeomtric p-values for overrepresentation of each GO term
 in the specified category among the GO annotations for the interesting
 genes (as indicated by the Entrez Gene ids).
}
\usage{
GOHyperG(x, lib, what="MF", universe=NULL)
}
\arguments{
  \item{x}{A character vector of unique Entrez Gene identifiers. }
  \item{lib}{The name of the annotation data package for the chip that
    was used or \code{"YEAST"}, see details for more information.}
  \item{what}{One of "MF", "BP", or "CC" indicating which of the GO
    categories to use for the computation.  In \code{GOKEGGHyperG},
    what can also be "KEGG"}
  \item{universe}{A character vector of unique Entrez Gene identifiers
    or \code{NULL}.  This is the population (the urn) of the
    Hypergeometric test.  When \code{NULL} (default), the population is
    all Entrez Gene ids in the annotation package that have a GO term
    annotation in the specified GO category (see details).}
}
\details{
  The Entrez Gene ids given in \code{x} define the selected set of
  genes.  The universe of Entrez Gene ids is determined by the chip
  annotation data package (\code{lib}) or specified by the
  \code{universe} argument which must be a subset of the Entrez Gene ids
  represented on the chip.  Both the selected genes and the universe are
  reduced by removing Entrez Gene ids that do not have any annotations
  in the specified GO category.

  For each GO term in the specified category that has at least one
  annotation in the selected gene set (\code{x}), we determine how many
  of its Entrez Gene annotations are in the universe set and how many
  are in the selected set.  With these counts we perform a
  Hypergeometric test using \code{phyper}.  This is equivalent to using
  Fisher's exact test.

  It is important that the correct chip annotation data package be
  identified as it determines the GO term to Entrez Gene id mapping as
  well as the universe of Entrez Gene ids in the case that the
  \code{universe} argument is omitted.

  For S. cerevisiae if the \code{lib} argument is set to \code{"YEAST"}
  then comparisons and statistics are computed using common names
  and are with respect to all genes annotated in the S. cerevisiae genome
  not with respect to any microarray chip.  This will \bold{not} be the
  right thing to do if you are working with a yeast microarray.
}

\note{
  Typically, one has a set of interesting genes/probes obtained from a
  microarray experiment and is interested in determining whether there
  is an overrepresentation of these genes at particular GO terms.
  \code{GOHyperG} carries out simple Hypergeometric tests to assess the
  overrepresentation of GO terms.

  
  Two substantial issues arise.  First, it is not clear how to do any
  form of p-value correction.  The tests are not independent and the
  underlying structure of the GO graph presents certain problems that
  need to be addressed.  The second substantial issue is that not all
  probes on a microarray map to a unique Entrez Gene identifer.  In
  \code{GOHyperG} every attempt to appropriately correct for
  non-uniqueness of mappings has been made.
}

\value{
  The returned value is a list with components:
  \item{pvalues }{The ordered p-values.}
  \item{goCounts}{The vector of counts of Entrez Gene ids from the
    universe at each node.}
  \item{intCounts}{The vector of counts of the supplied Entrez Gene ids
    annotated at each GO term.}
  \item{numLL}{The number of unique Entrez Gene ids in the universe that
    are mapped to some term in the specified GO category.}
  \item{numInt}{The number of unique Entrez Gene ids in the selected
    gene set, \code{x}, that are mapped to some term in the specified GO
    category.}
  \item{chip}{A string identifying the chip annotation data package used.}
  \item{intLLs}{The input vector \code{x}.}
  \item{go2Affy}{A list with one element for each GO term tested, containing
    the Affymetrix identifiers associated with that node, for the whole
    chip (not just the interesting genes).  This is the same as
    extracting the tested GO ids from the annotation package's
    GO2ALLPROBES environment.}
}
\author{R. Gentleman}

\seealso{
  \code{\link{hyperGTest}}
  \code{\link[Category]{geneKeggHyperGeoTest}}
  \code{\link[Category]{geneCategoryHyperGeoTest}}
  \code{\link{phyper}}
}

\examples{
\dontrun{
library(hgu95av2.db)
library(GO)
w1<-as.list(hgu95av2ENTREZID)
w2<-unique(unlist(w1))
set.seed(123)
## pick a 25 interesting genes
myLL <- sample(w2, 25)
xx <- GOHyperG(myLL, lib="hgu95av2.db", what="CC")
xx$numLL
xx$numInt
sum(xx$pvalues < 0.01)
}
}
\keyword{ htest }
