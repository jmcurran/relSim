\name{readFreqs}
\alias{readFreqs}
\title{Read in a file of allele frequencies}
\description{Reads in a file of alleles in a particular format.
}
\usage{
readFreqs(strPath)
}
\arguments{
  \item{strPath}{The file from which to read the frequencies}
}
\value{
  a list containing two lists, loci and freqs. The list loci is a vector
  of the locus names in the frequency file. The list freqs, is a list of
  vectors with each vector containing the frequencies of the alleles at
  the locus. The names of the elements of the vectors are the STR allele designations.
}
\details{
  This function reads frequencies in a very particular text format. The
  first line contains the number of loci in the multiplex. The next line
  will contain the name of the first locus and the number of alleles, nA,
  the locus separated by a comma. The next nA lines contain the allele
  number (from 1 to nA), the STR designation of the allele, and the
  frequency separated by commas. This pattern is repeated for each
  locus. In the future this function will read the rectangular allele
  freqency table used by FSI Genetics and other journals.
}
\author{James M. Curran}
