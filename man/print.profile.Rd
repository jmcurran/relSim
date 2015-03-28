\name{print.profile}
\alias{print.profile}
\title{Print a DNA profile}
\description{
  Nicely prints a profile object out in genotype pairs
}
\usage{
\method{print}{profile}(x, horizontal = FALSE, \dots)
}
\arguments{
  \item{x}{The profile object to be printed}
  \item{horizontal}{if \code{TRUE} then the profile will print on a
    single line instead of multiple lines. Useful for comparing two
    profiles}
  \item{\dots}{Ignored - really should be passed to print, but given cat
    is actually called they are ignored}
}
\examples{
data(fbiCaucs)
P1 = randomProfile(fbiCaucs)
P2 = randomProfile(fbiCaucs)
P1
print(P1, horizontal = TRUE)
print(P2, horizontal = TRUE)
}
\author{James M. Curran}
