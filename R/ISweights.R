ISW = function(){
  data("USCaucs")
  f = USCaucs$freqs[[1]]
  x = IS(f, 100, 4, 2)
  X = t(apply(x, 1, function(row)which(row > 0)))
  
  library(multicool)
  m = initMC(c(1,2))
  perms = allPerm(m)
  
  
  p = ISprobs(f, X, perms)
  return(p)
}