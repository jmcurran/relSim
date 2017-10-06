#' @importFrom multicool initMC allPerm
ISW = function(N = 1e5, numContrib = 4, numPeaks = 3){
  data("USCaucs")
  f = USCaucs$freqs[[2]]
  r = IS(f, N, numContrib, numPeaks)
  x = r$Alleles
  X = apply(x, 1, function(row){
    alleles = which(row > 0)
    counts = row[alleles]
    freqs = f[alleles]
    return(list(a = alleles, c = counts, f = freqs))
  })
  
  if(numPeaks != 1){
    m = initMC(1:numPeaks)
    perms = allPerm(m)
  }else{
    perms = matrix(1, 1, 1)
  }
  
  w = ISprob(X, perms)
  p = mean(exp(r$probs - w))
  
  return(list(Alleles = X, perms = perms, weights = w, probs = r$probs, p = p))
}