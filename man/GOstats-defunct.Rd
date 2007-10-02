\name{GOstats-defunct}
\alias{GOstats-defunct}
\title{Defunct Functions in GOstats Package}
%------ PLEASE: put \alias{.} here for  EACH !
\alias{combGOGraph}
\alias{hyperGtable}
\alias{hyperG2Affy}
\description{
  The functions or variables listed here are no longer part of GOstats as
  they are not needed (any more).
}
\usage{
combGOGraph(g1, g2)
hyperGtable(probids, lib, type="MF", pvalue=0.05,
                        min.count=10, save = FALSE, output = TRUE,
                        filename = NULL, universe = NULL)
hyperG2Affy(probids, lib, type="MF", pvalue=0.05,
                         min.count=10, universe = NULL)

}
\arguments{
  \item{g1}{A graph}
  \item{g2}{A graph}
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
  \item{save}{Boolean - Set to \code{TRUE} to save the resulting \code{data.frame}.}
  \item{output}{ Boolean - Set to \code{TRUE} to output the resulting
    \code{data.frame} as a text file.}
  \item{filename}{ If output is set to \code{TRUE}, give the file name
    for the output file.}
  
  }
\details{
  \code{combGOGraph} was replaced by \code{join}.
  \code{hyperGtable} was replaced by \code{summary}.
  \code{hyperG2Affy} was replaced by \code{probeSetSummary}.
}

\keyword{manip}
