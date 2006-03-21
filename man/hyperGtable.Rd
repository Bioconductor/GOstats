\name{hyperGtable}
\alias{hyperGtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Table from GOHyperG }
\description{
  Output a table from a call to \code{\link[GOstats]{GOHyperG}}.
}
\usage{
hyperGtable(probids, lib, type = "MF", pvalue = 0.05, min.count = 10,
save = FALSE, output = TRUE, filename = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probids}{ A vector of Affymetrix probe IDs}
  \item{lib}{ An annotation package (e.g., \code{hgu95av2})}
  \item{type}{One of "MF", "CC", "BP", indicating molecular function,
    cellular component, or biological process, respectively.}
  \item{pvalue}{The significance level used to choose GO terms}
  \item{min.count}{ The minimum number of a given GO term that must be
    on the chip in order to choose that GO term. This protects against
    very low p-values that result from the situation where there are
    very few genes with a given GO term on the chip, but one or two are
    found in the set of significant genes.}
  \item{save}{Boolean - Set to \code{TRUE} to save the resulting \code{data.frame}.}
  \item{output}{ Boolean - Set to \code{TRUE} to output the resulting
    \code{data.frame} as a text file.}
  \item{filename}{ If output is set to \code{TRUE}, give the file name
    for the output file.}
}
\details{
 This function will take a vector of Affymetrix IDs and output a
 \code{data.frame} containing information on all GO terms that are found
 significant in a call to \code{\link[GOstats]{GOHyperG}}.
}
\value{
  A \code{data.frame} containing the following columns:

  GO Term - The GOID associated with a given GO term
  
  Description - The description for that GO term
  
  p-values - The p-value associated with that GO term
  
  Number significant - The number of genes in the vector of Affymetrix IDs
  with that particular GO term
  
  Number on chip - The number of genes on the Affymetrix chip with that
  particular GO term.
}
\author{ James W. MacDonald <jmacdon@med.umich.edu>}
\examples{
library("hgu95av2")
probids <- ls(hgu95av2GENENAME)[1:10]
hyperGtable(probids, "hgu95av2", save = TRUE, output = FALSE)
}
\keyword{ manip }% __ONLY ONE__ keyword per line
\keyword{htest}

