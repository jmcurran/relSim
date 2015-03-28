\name{randomProfile}
\alias{randomProfile}
\title{Generate a random DNA profile from a given set of allele frequencies}
\description{
  Generates a random DNA profile from a given set of allele frequencies.
}
\usage{
randomProfile(Freqs)
}
\arguments{
  \item{Freqs}{A list containing two lists labelled loci and freqs. The
second list is a list of vectors containing the allele frequencies of
each allele at each locus in the multiplex.} 
}
\value{
  A matrix with 2 columns and nLoci rows. Each row represents the
  genotpe of the random individual at that locus. The genotype alleles
  are always ordered so that allele1 <= allele2.
}
\details{
  The alleles are simply integers rather than the STR repeat
  numbers. This speeds up computation immensely when calculating any of
  the LRs or IBS.
}
\seealso{
  randomChild, randomSample, randomSib
}
\examples{
data(fbiCaucs)
P1 = randomProfile(fbiCaucs)
}
\author{James M. Curran}
