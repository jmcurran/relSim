#' Read a set of profiles from a file
#'
#' Reads a set of profiles from a file
#'
#' The alleles are recorded integers rather than the STR repeat numbers. This
#' speeds up computation immensely when calculating any of the LRs or IBS.
#'
#' @param fileName a path to the profile file.
#' @param Freqs A list containing two lists labelled loci and freqs. The second
#' list is a list of vectors containing the allele frequencies of each allele
#' at each locus in the multiplex.
#' @param delim a character that delimits the fields in the profile file.
#' @param header a boolean witch is \code{TRUE} if the profile file has a column header line.
#' @return matrix of profiles.
#' @author James M. Curran
#' @examples
#'
#' data(fbiCaucs)
#'
#' @export readProfiles
readProfiles = function(fileName,
                        freqs,
                        delim = "\t",
                        header = FALSE,
                        id = 1) {
  if (header) {
    prof = read.table(fileName,
                      sep = delim,
                      header = FALSE,
                      skip = 1)
  } else{
    prof = read.table(fileName, sep = delim, header = FALSE)
  }
  
  idCol = NULL
  if (id != -1) {
    idCol = prof[, id]
    prof = prof[, -id]
  }
  numLoci = length(freqs$freqs)
  
  for (loc in 1:numLoci) {
    locus = freqs[["freqs"]][[loc]]
    Alleles = as.numeric(names(locus))
    
    i1 = 2 * loc - 1
    i2 = i1 + 1
    prof[, i1] = match(prof[, i1], Alleles)
    prof[, i2] = match(prof[, i2], Alleles)
  }
  
  rownames(prof) = 1:nrow(prof)
  
  return(prof)
}
