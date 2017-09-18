ISW = function(N = 1e5, numContrib = 4, numPeaks = 3){
  data("USCaucs")
  f = USCaucs$freqs[[1]]
  r = IS(f, N, numContrib, numPeaks)
  x = r$Alleles
  X = t(apply(x, 1, function(row)which(row > 0)))
  
  m = initMC(1:numPeaks)
  perms = allPerm(m)
  
  w = ISprob(f, X, perms)
  p = exp(r$probs) / w
  
  return(list(Alleles = r$Alleles, perms = perms, weights = w, probs = exp(r$probs)))
}