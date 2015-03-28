\name{toNexus}
\alias{toNexus}
\title{Export a population with substructure to a Nexus file}
\description{Exports a population with population substructure to a
  Nexus formatted file so that GDA can be used to check the Fst calculations
}
\usage{
toNexus(Pop, fileName = 'output.nex')
}
\arguments{
  \item{Pop}{An object of type 'population' - see \code{breedFst} for a
    description of the object}
  \item{fileName}{The name of the file output file}

}
\references{
  Maddison DR, Swofford DL, Maddison WP (1997), NEXUS: An extensible
  file format for systematic information, Systematic Biology 46 (4):
  590--621.

  Zaykin, D. and Lewis, P., GDA - software to accompany Genetic Data
  Analysis II, \url{http://hydrodictyon.eeb.uconn.edu/people/plewis/software.php}.
}
\examples{
data(USCaucs)
p = breedFst(USCaucs)
toNexus(p)
}
\seealso{breedFst}
\author{James M. Curran}

