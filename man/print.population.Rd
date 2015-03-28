\name{print.population}
\alias{print.population}
\title{Print summary details of a substructed population}
\description{
  Nicely prints summary information about a substructured population
  created using \code{breedFst}
}
\usage{
\method{print}{population}(x, \dots)
}
\arguments{
  \item{x}{The population object to be printed}
  \item{\dots}{Ignored - really should be passed to print, but given cat
    is actually called they are ignored}
}
\examples{
data(fbiCaucs)
p = breedFst(fbiCaucs)
print(p)
}
\seealso{breedFst}
\author{James M. Curran}
