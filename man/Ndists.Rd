\name{Ndists}
\alias{Ndists}
\alias{Bdists}
\docType{data}
\title{Distance matrices for the BCR/ABL and NEG subgroups.}
\description{
  These are precomputed distance matrices between all transcription
  factors selected. In the future they will be computed on the fly but
  currently that takes about 3 hours and so precomputed versions are
  supplied. 
}
\usage{data(Ndists)
  data(Bdists)}
\format{
  These are both distance matrices.
}
\source{
 They are based on the ALL data, \code{\link{ALL}}.
}
\examples{
data(Ndists)
data(Bdists)
}
\keyword{datasets}
