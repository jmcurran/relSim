\name{checkFreqs}
\alias{checkFreqs}
\title{Make sure that the frequencies are such}
\description{
  Checks whether a list of frequencies at a series of genetic loci both
  sum to one and lie between 0 and 1.
}
\usage{
checkFreqs(Freqs)
}
\arguments{
  \item{Freqs}{A list containg elements \code{loci} and
    \code{freqs}. \code{freqs} is a list of vectors containing the
    frequencies at the given loci.}
}
\details{
  If a locus fails to sum to one, or there are alleles which fall below
  zero or above one, then a warning message will be returned for each
  item in error.
}
\seealso{
  normalizeFreqs
}
\examples{
data(fbiCaucs)
checkFreqs(fbiCaucs)

## induce an error
fbiCaucs$freqs[[1]] = runif(10)
checkFreqs(fbiCaucs)
}
\author{James M. Curran}
