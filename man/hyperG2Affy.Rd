\name{hyperG2Affy}
\alias{hyperG2Affy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{(DEPRECATED) A List of Affymetrix IDs from GOIDs}
\description{
  Use \code{probeSetSummary} instead.
  
  Output a list of Affymetrix IDs associated with each significant GOID,
  based on the output from hyperGtable.
}
\usage{
hyperG2Affy(probids, lib, type = "MF", pvalue = 0.05, min.count = 10,
universe = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probids}{ A vector of Affymetrix probe IDs}
  \item{lib}{ An annotation package (e.g., \code{hgu95av2.db})}
  \item{type}{One of "MF", "CC", "BP", indicating molecular function,
    cellular component, or biological process, respectively.}
  \item{pvalue}{The significance level used to choose GO terms}
  \item{min.count}{ The minimum number of a given GO term that must be
    on the chip in order to choose that GO term. This protects against
    very low p-values that result from the situation where there are
    very few genes with a given GO term on the chip, but one or two are
    found in the set of significant genes.}
   \item{universe}{A character vector of unique Entrez Gene identifiers
    or \code{NULL}.  This is the population (the urn) of the
    Hypergeometric test.  When \code{NULL} (default), the population is
    all Entrez Gene ids in the annotation package that have a GO term
    annotation in the specified GO category (see \code{GOHyperG} for
    more details).}
}
\details{
 The function \code{hyperGtable} will output a table listing a set of
 GOIDs that are significantly over-represented in a set of significant
 genes. The next logical step might be to look at the individual genes
 that are represented by each GOID. This function will output the
 Affymetrix IDs for each GOID in a list, which can then be used to
 either output HTML tables using \code{probes2table}, or possibly just
 output the gene names for each GOID. 
}
\value{
  A named \code{list} with all the Affymetrix IDs associated with each
  significant GOID.
}
\author{ James W. MacDonald <jmacdon@med.umich.edu>}
\examples{
\dontrun{
library("hgu95av2.db")
probids <- ls(hgu95av2GENENAME)[1:10]
hyperG2Affy(probids, "hgu95av2.db")
}
}
\keyword{ manip }% __ONLY ONE__ keyword per line
\keyword{htest}

