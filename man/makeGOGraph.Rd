\name{makeGOGraph}
\alias{makeGOGraph}
\title{ Construct a GO Graph }
\description{
  The directed acyclic graph (DAG) based on finding the most specific
  terms for the supplied Entrez Gene IDs is constructed and returned. The
  constructuion is per GO ontology (there are three, MF, BP and CC) and
  once the most specific terms have been identified then all less
  specific terms are found (these are the parents of the terms) and then
  their parents and so on, until the root is encountered.
}
\usage{
makeGOGraph(x, Ontology = "MF", removeRoot = TRUE, mapfun = NULL,
            chip = NULL)
}
\arguments{
  \item{x}{A vector of Entrez Gene IDs. }
  \item{Ontology}{Which GO ontology to use (CC, BP, or MF). }
  \item{removeRoot}{A logical value indicating whether the GO root node
    should be removed or not.}

  \item{mapfun}{A function taking a character vector of Entrez Gene IDs
    as its only argument and returning a list of "GO lists" matching the
    structure of the lists in the GO maps of annotation data packages.
    The function should behave similarly to
    \code{mget(x, eg2gomap, ifnotfound=NA)},
    that is, \code{NA} should be returned if a specified Entrez ID has
    no GO mapping.  See details for the interaction of \code{mapfun} and
    \code{chip}.}

  \item{chip}{The name of a DB-based annotation data package (the name
    will end in ".db").  This package will be used to generate an Entrez
    ID to GO ID mapping instead of \code{mapfun}.}

}

\details{
  For each supplied Entrez Gene identifier all the GO annotations (in the
  specified ontology) are found.  The mapping is achieved in one of
  three ways:

  \enumerate{
    \item If \code{mapfun} is provided, it will be used to perform the
      needed lookups.  In this case, \code{chip} will be ignored.

    \item If \code{chip} is provided and \code{mapfun=NULL}, then the
      needed lookups will be done based on the Entrez to GO mappings
      encapsulated in the specified annotation data package.  This is
      the recommended usage.

    \item If \code{mapfun} and \code{chip} are \code{NULL} or missing,
      then the function will attempt to load the GO package (the
      environment-based package, distinct from GO.db).  This package
      contains a legacy environment mapping Entrez IDs to GO IDs.  If
      the GO package is not available, an error will be raised.
      Omitting both \code{mapfun} and \code{chip} is not recommended as
      it is not compatible with the DB-based annotation data packages.
  }

  The mappings are different for the different ontologies. Typically a
  GO indentifier is used only in one specific ontology.

  The resulting structure is stored in a graph using the \code{graph}
  package, again from Bioconductor.

}
\value{
  An object that inherits from the \code{graph} class. The particular
  implementation is not specified.
}
\references{The Gene Ontology Consortium }
\author{R. Gentleman}

\seealso{\code{\link{oneGOGraph}}}

\examples{
 library("hgu95av2.db")
 set.seed(321)
 gN <- unique(sample(keys(hgu95av2.db, 'ENTREZID'), 4))
 gg1 <- makeGOGraph(gN, "BP", chip="hgu95av2.db")

}
\keyword{manip}
