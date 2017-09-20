#' @importFrom multicool initMC allPerm
ISW = function(N = 1e5, numContrib = 4, numPeaks = 3){
  data("USCaucs")
  f = USCaucs$freqs[[1]]
  r = IS(f, N, numContrib, numPeaks)
  x = r$Alleles
  X = apply(x, 1, function(row){
    alleles = which(row > 0)
    counts = row[alleles]
    freqs = f[alleles]
    return(list(a = alleles, c = counts, f = freqs))
  })
  newFreqs = t(apply(x, 2, function(row)f[which(row > 0)]))
  newFreqs = newFreqs / rowSums(newFreqs)
  
  m = initMC(1:numPeaks)
  perms = allPerm(m)
  
  w = ISprob(f, X, newFreqs, perms)
  p = exp(r$probs) / w
  
  return(list(Alleles = r$Alleles, freqs = newFreqs, perms = perms, weights = w, probs = exp(r$probs)))
}