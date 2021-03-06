% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lrPC.R
\name{lrPC}
\alias{lrPC}
\title{Likelihood Ratio for Parent-Child / Paternity Index}
\usage{
lrPC(parent, child, Freqs = NULL, nLoci = length(parent)/2, f = NULL, n = NULL)
}
\arguments{
\item{parent}{A matrix consisting of 2 columns and nLoci rows. Each entry in
the matrix is the (coded) allele held by the individual. This represents the
alleged parent. The relationship is reflexive so it does not matter which
profile is labelled parent and child.}

\item{child}{See \code{parent}}

\item{Freqs}{A list containing two lists labelled loci and freqs. The second
list is a list of vectors containing the allele frequencies of each allele
at each locus in the multiplex. This argument or both f and n must be
specified}

\item{nLoci}{The number of loci in the profiles}

\item{f}{A concatenated vector of allele frequencies. Specifying this speeds
up computation enormously}

\item{n}{A vector of length \code{nLoci} giving the number of alleles at
each locus. Specifying this in advance enormously speeds up computation}
}
\value{
A value between 0 and infinity representing support (or lack of
support if the value is less than 1) for the hypothesis that the two
profiles are parent and child. There is no mutation built into this
calculation. This means that the LR will be zero if the profiles do not
share at least one allele in common at each locus in the multiplex.
}
\description{
Calculates Likelihood Ratio comparing the probability of two profiles if
they are indeed parent-child compared to unrelated. This is the paternity
index or PI.
}
\examples{

data(fbiCaucs)
P1 = randomProfile(fbiCaucs)
C1 = randomChild(P1, fbiCaucs)
lrPC(P1, C1, fbiCaucs)

}
\references{
Buckleton, J, Triggs, C.M., and Walsh, S.J. (2005)\emph{Forensic
DNA Evidence Interpretation}, CRC Press., Boca Raton, FL. p.410
}
\seealso{
lrSib, IBS
}
\author{
James M. Curran
}
