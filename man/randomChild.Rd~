\name{randomChild}
\alias{randomChild}
\title{Generate a random child from a given DNA profile and
  a given set of allele frequencies}
\description{
  Generates a random child (or parent) from a given DNA profile from a
  given set of allele frequencies. At each locus, the child
  inherits the first allele of the given profile with one half, or the
  second allele with probability one half. The second allele is chosen
  at random with probability proportional to the allele frequencies.
}
\usage{
randomChild(profile, Freqs)
}
\arguments{
  \item{profile}{A matrix consisting of 2 columns and nLoci rows. Each entry
in the matrix is the (coded) allele held by the individual. This
represents the parent. The relationship is reflexive so it
does not matter if the profile is a parent or a child.}
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
C1 = randomChild(P1,fbiCaucs)
P1
C1
}
\author{James M. Curran}
