\name{randomSibPairs}
\alias{randomSibPairs}
\title{Generate one or more pairs of random siblings from a given set of
  allele frequencies
}
\description{
  Generates one or more pairs of random siblings from a given set of
  allele frequencies.
}
\usage{
randomSibPairs(Freqs, BlockSize = 1)
}
\arguments{
  \item{Freqs}{A list containing two lists labelled loci and freqs. The
second list is a list of vectors containing the allele frequencies of
each allele at each locus in the multiplex.}
\item{BlockSize}{The number of pairs of profiles to generate}
}
\value{
  A list of length \code{BlockSize}. Each element of the list has a
  sublist containing two profiles called \code{sib1} and \code{sib2}
}
\details{
  The alleles are simply integers rather than the STR repeat
  numbers. This speeds up computation immensely when calculating any of
  the LRs or IBS.
}
\seealso{
  randomPCPairs, randomProfilePairs
}
\examples{
data(fbiCaucs)
P = randomSibPairs(fbiCaucs)
P$sib1
P$sib2
}
\author{James M. Curran}
