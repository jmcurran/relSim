\name{fbiCaucs}
\docType{data}
\alias{fbiCaucs}
\title{CODIS STR Loci allele frequency data}
\description{
  This data structure
}
\usage{
data(fbiCaucs)
}
\format{
  This data set is a list which has two sub-lists. The lists are named
  loci and freqs. loci is a vector of the 13 CODIS STR locus
  names. freqs is a list of 13 vectors, each vector contains the allele
  frequencies published for US Caucasians in Budowle et al. (2001).
}
\references{
  Budowle B, Shea B, Niezgoda S, Chakraborty R. (2001), \emph{CODIS STR
    loci data from 41 sample populations}, J. Forensic Sci. 46:453-89.
}
\seealso{
  USCaucs
}
\examples{
data(fbiCaucs)
names(fbiCaucs)
fbiCaucs$loci
names(fbiCaucs$freqs)
fbiCaucs$freqs[[1]]
names(fbiCaucs$freqs[[1]])
fbiCaucs$freqs[[1]][1]
}
\author{James M. Curran}
\keyword{datasets}
