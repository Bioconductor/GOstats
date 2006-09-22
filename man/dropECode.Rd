\name{dropECode}
\alias{dropECode}
\title{ Drop GO labels for specified Evidence Codes }
\description{
  Genes are mapped to GO terms on the basis of evidence codes. In some
  analyses it will be appropriate to drop certain sets of annotations
  based on specific evidence codes.
}
\usage{
dropECode(inlist, code="IEA")
}
\arguments{
  \item{inlist}{A list of GO data }
  \item{code}{The set of codes that should be dropped. }
}
\details{
  A simple use of \code{\link{lapply}} and \code{\link{sapply}} to find
  and eliminate those terms that have the specified evidence codes.

  This might be used when one is using to GO to validate a sequence
  matching experiment (for example), then all terms whose mapping was
  based on sequence similarity (say ISS and IEA) should be removed.
  
}
\value{
  A list of the same length as the input list retaining only those
  annotations whose evidence codes were not the ones in the exclusion
  set \code{code}.
}
\author{R. Gentleman }

\seealso{\code{\link{getEvidence}}, \code{\link{getOntology}} }

\examples{
 bb=GOLOCUSID2GO$"4121"
 getEvidence(bb[1:3])
 cc = dropECode(bb[1:3])
 if (length(cc))
   getEvidence(cc)
}


}
\keyword{manip }
