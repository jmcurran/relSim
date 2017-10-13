#' @importFrom multicool initMC allPerm
#' @export
IS = function(freqs, numContributors = 4, numIterations = 100){
  if(numContributors < 2){
    stop("numContributors must be >= 2")
  }
  
  m = initMC(1:4)
  perms = list(allPerm(m))
  r = .IS(USCaucs$freqs, numIterations, 4, 4, perms)
  # r = .IS(f, N, numContrib)
  # x = r$Alleles
  # X = apply(x, 1, function(row){
  #   alleles = which(row > 0)
  #   counts = row[alleles]
  #   freqs = f[alleles]
  #   n = length(alleles)
  #   return(list(a = alleles, c = counts, f = freqs, n = n))
  # })
  # 
  # perms = vector(length = maxPeaks, mode = "list")
  # for(np in 1:maxPeaks)
  # if(np != 1){
  #   m = initMC(1:np)
  #   perms[[np]] = allPerm(m)
  # }else{
  #   perms[[np]] = matrix(1, 1, 1)
  # }
  # 
  # w = ISprob(X, perms) - log(maxPeaks)
  # p = mean(exp(r$probs - w))
  # 
  # return(list(Alleles = X, perms = perms, weights = w, probs = r$probs, est = p))
  m = initMC(1:np)
  return(r)
}

ISprobR = function(X, perms){
  freqs = X$f
  numAlleles = X$n
  numPerms = nrow(perms)
  
  r = 0
  
  for(j in 1:numPerms){
    p = s = freqs[perms[j, 1]];
    
    for(k in 2:numAlleles){
      pk = freqs[perms[j, k]];
      p  = p *  pk / (1 - s);
      s  = s + pk;
    }
    
    r = r + p
  }
  
  r = log(r)
  
  freqs = freqs / sum(freqs)
  counts = X$c - 1
  p2 = sum(counts * log(freqs) - log(factorial(counts))) + log(factorial(sum(counts)))
  
  return(r + p2) 
}