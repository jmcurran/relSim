\name{fetchBMdata}
\alias{fetchBMdata}
\title{Retrieve data from Budowle and Moretti (1999) from the web}
\description{
  Retreives the Budowle and Moretti (1999) and compiles the allele
  frequency tables needed for the other parts of this package such as \code{sim}.
}
\usage{
fetchBMdata()
}
\value{
  A list consisting of six elements corresponding to the six populations
  detailed in the data set. Each of the list elements is a list in
  itself with two further elements named loci and freqs. loci is a
  vector of the 13-20 STR locus names. freqs is a list of 13-20 vectors,
  each vector contains the allele frequencies.
}
\details{
  The first three populations have data on 20 loci, the second three on
  13 loci. The missing values (0's in the raw data) have been dropped
  and are not used in calculating the frequencies. This function will
  not work if you are not connected to the internet, or access to the
  internet is blocked.
}
\references{
  Budowle, B. and Moretti, T.R. (1999), \emph{Genotype Profiles for Six
    Population Groups at the 13 CODIS Short Tandem Repeat Core Loci and
    Other PCR Based Loci}, Forensic Science Communications 1(2).
}
\seealso{
  fbiCaucs, USCaucs
}
\examples{
## not run
\dontrun{
db = fetchBMdata()
names(db)
f = db[["TRINIDADIAN"]]$freqs
dbExpect(f, k = "UN", collapse = TRUE)
}
}
\author{James M. Curran}
\keyword{datasets}
