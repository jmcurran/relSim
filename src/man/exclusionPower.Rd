% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exclusionPower.R
\name{exclusionPower}
\alias{exclusionPower}
\alias{ep}
\title{Calculate the exclusion power of a multiplex by locus}
\usage{
exclusionPower(Freqs)
}
\arguments{
\item{Freqs}{A list containing two vectors and a list, called loci, counts,
and freqs. The elements of loci are the loci present in the multiplex. The
elements are freqs a vectors of allele frequencies for the locus. The
elements of counts are irrelevant here.}
}
\value{
The exclusion power for each locus.
}
\description{
Calculates the exclusion power
\deqn{1-2\left(\sum_{i=1}^{n_l}p_i^2\right)^2-4\sum_{i=1}^{n_l}p_i^4}{1-2*sum(pi^2)^2-4*sum(pi^4)}
at each locus for a set of allele frequencies.
}
\examples{

data(USCaucs)
ep(USCaucs)

## get the multiplex wide exclusion power
1 - prod(1-ep(USCaucs))

}
\references{
NRC II, Evaluation of Forensic Evidence, (1996), p.96, National
Academy Press.
}
\author{
James M. Curran
}
