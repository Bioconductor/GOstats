\name{idx2dimnames}
\alias{idx2dimnames}
\title{Index to Dimnames}
\description{
  A function to map from integer offsets in an array to the
corresponding values of the row and column names. There is probably a
better way but I didn't find it.
}
\usage{
idx2dimnames(x, idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{matrix} or \code{data.frame}. }
  \item{idx}{An integer vector of offsets into the matrix (values
  between 1 and the \code{length} of the matrix.}
}

\value{
 A list with two components.
  If it is a LIST, use
  \item{rowNames }{The row names corresponding to the integer index.}
  \item{colNames }{The column names corresponding to the integer index.}
}

\author{R. Gentleman}

\seealso{\code{\link{dimnames}} }
\examples{
 data(Ndists)
 ltInf = is.finite(Ndists)
 xx = idx2dimnames(Ndists, ltInf)
}
\keyword{manip}
