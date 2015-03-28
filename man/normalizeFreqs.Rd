\name{normalizeFreqs}
\alias{normalizeFreqs}
\title{Normalize frequencies to 1}
\description{
  Normalize a list of frequencies at a series of genetic loci both
  sum to one. Not that this does not deal with the problem of values
  larger than one or smaller than zero.
}
\usage{
normalizeFreqs(Freqs)
}
\arguments{
  \item{Freqs}{A list containg elements \code{loci} and
    \code{freqs}. \code{freqs} is a list of vectors containing the
    frequencies at the given loci.}
}
\details{
  Divides vector in Freqs$freqs by the vector sum.
}
\value{
  A list containg elements \code{loci} and
    \code{freqs}. \code{freqs} is a list of vectors containing the
    frequencies at the given loci.
}
\seealso{
  checkFreqs
}
\examples{
data(fbiCaucs)

## induce an error
fbiCaucs$freqs[[1]] = rgamma(10,1,1)
checkFreqs(fbiCaucs)
fbiCaucs = normalizeFreqs(fbiCaucs)
checkFreqs(fbiCaucs)
}
\author{James M. Curran}
