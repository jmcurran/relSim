rdirichlet = function(N, counts){
  nRow = length(counts)
  nCol = N
  X = rgamma(nRow * nCol, shape = rep(counts, rep(nCol, nRow)), rate = 1)
  X = matrix(X, nrow = nRow, byrow = TRUE)
  X = prop.table(X, 2)
  return(X)
}

cv = function(x){sd(x)/mean(x)}
neff = function(x){N / (1 + cv(x)^2)}

#library(relSim)
#Freqs = readFreqs("D:/Dropbox/Dropbox/Work/2015/06-June/Coble, Mike/Data/alln1.csv")
SE33 = c(5e-04,5e-04,5e-04,5e-04,0.001,0.0043,0.0019,0.0116,0.0019,0.0319,0.0034,0.0376,0.0024,0.0487,0.0019,0.001,0.0734,5e-04,0.0019,0.096,0.0931,0.0029,0.0724,0.0068,0.0333,0.0179,0.0145,0.0232,0.0043,0.028,5e-04,0.0232,0.0352,5e-04,0.0589,5e-04,0.0738,5e-04,0.0642,0.001,0.0468,5e-04,0.0376,0.0014,0.0188,5e-04,0.0092,0.001,0.0034,0.0039,5e-04,5e-04)
names(SE33) = c(6.3,7.0,10.2,11.0,11.2,12.0,12.2,13.0,13.2,14.0,14.2,15.0,15.2,16.0,16.2,16.3,17.0,17.2,17.3,18.0,19.0,19.2,20.0,20.2,21.0,21.2,22.0,22.2,23.0,23.2,24.0,24.2,25.2,26.0,26.2,27.0,27.2,27.3,28.2,28.3,29.2,30.0,30.2,31.0,31.2,32.0,32.2,33.0,33.2,34.0,34.2,36.0)
#SE33 = Freqs$freqs$SE33
log.SE33 = log(SE33)
SE33.counts = floor(SE33 * 2072)
nA = length(SE33)


#set.seed(5284771)
modifyProbs = function(probs, nBig = 6, pTarget = 0.75, plot = FALSE){
  o = rev(order(probs))
  big = o[1:nBig]
  
  ## the nBiggest alleles get pTarget of the mass, the rest get 1-pTarget 
  u = rgamma(nBig, 50, 1)
  u = sort(u / sum(u), decreasing = TRUE)
  modProb = probs
  
  i = 1
  for(b in big){
    modProb[b] = u[i] * pTarget
    i = i + 1
  }
  
  scaleFactor = (1 - pTarget) / sum(modProb[-big]) 
  modProb[-big] = modProb[-big]  * scaleFactor
  if(plot)
    barplot(rbind(probs, modProb), beside = T)
  return(modProb)
}

symmStat = function(x, q = 0.05){
  qtls = quantile(x, probs = c(q * 0.5, 0.5, 1 - q * 0.5))
  return((qtls[2] - qtls[1]) / (qtls[3] - qtls[1]))
}

weightSim = function(){
  P = seq(0.50, 0.90, length = 41)
  k = 1
  s = m = rep(0, length(P))
  pb = txtProgressBar(0, length(P), 0, style = 3)
  for(p in P){
    newProbs = modifyProbs(SE33, 6, p)
    w = relSim:::.sampleWeights(SE33, newProbs, 6, 1e4)
    n = relSim:::.tabulateN(newProbs, 6, 1e6)
    s[k] = symmStat(w)
    m[k] = sum(n * 1:12) / 1e6
    setTxtProgressBar(pb, k)
    k = k + 1
  }
  close(pb)
  m0 = sum((relSim:::.tabulateN(SE33, 6, 1e6) * 1:12) / 1e6)
  return(list(p = P, s = s, m = m, m0 = m0))
}

SE33.modProb = modifyProbs(SE33, 5)
SE33.modCounts = SE33.modProb * 2072


# X = rmultinom(1e7, 12, prob = SE33.modProb)
# w = exp(t(X) %*% matrix(log(SE33) - log(SE33.modProb),nc=1))
# nx = apply(X, 2, function(x)sum(x!=0))
# pt = c(100, 1000, 1e4, 1e5, seq(2e5, 1e6, length = 100))
# k = 1
# est = rep(0, length(pt))
# for(p in pt){
#   est[k] = mean(w[1:p] * (nx[1:p] <= 2))
#   k = k + 1
# }
# plot(pt, est, type = 'l')
# abline(h = 7.2e-9)
# 
# rm(X)
# table(nx)

sim = function(nRep = 100, N = 1e6, pTarget = 0.75, numContrib = 6, numBig = 6){
  log.SE33 = matrix(log(SE33), nr = 52, nc = chunkSize)
  
  pb = txtProgressBar(0, nRep, 0, style = 3)
  p1 = p2 = ne = rep(0, nRep)
  SE33.modProb = modifyProbs(SE33, numBig, pTarget)
  
  for(j in 1:nRep){
    p2[j] = relSim:::.importance(SE33, SE33.modProb, numContrib, N)
    setTxtProgressBar(pb, j)
    plot(p2[1:j], type = 'l')
    p3 = cumsum(p2[1:j]) / (1:j)
    lines(p3[1:j], col = 'red', lty = 3)
    
    if(any(p2[1:j] == 0)){
      idx = which(p2[1:j] == 0)
      points(idx, rep(0, length(idx)), pch = 'x', col = "red")
    }
    abline(h = 7.2e-9)
  }
  
  return(p2)
}

set.seed(123)
mcmcSim = function(N = 1000){
  g0 = rmultinom(1, 12, SE33)
  l0 = sum(g0 * log(SE33))
  log.SE33 = log(SE33)
  n0 = sum(g0 != 0)
  
  res = rep(0, N)
  toss = runif(N)
  log.u = log(runif(N))
  pb = txtProgressBar(0, N, 0, style = 3)
  cat(paste(l0, "\n"))
  for(i in 1:N){
    if(toss[i] < 0.5 && n0 < 12){
      ##cat(paste(i, n0, l0, " birth\n"))
      
      ## birth  - pick an allele not in genotype
      ## and decrease the count of one of the others that has count > 1
      g1 = g0
      SE33subs = SE33
      SE33subs[g1 != 0] = 0
      birth = sample(1:52, 1, prob = SE33subs)
      g1[birth] = 1
      k = which(g1 > 1)
      if(length(k) > 1)
        k = sample(k, size = 1)
      g1[k] = g1[k] - 1
      l1 = sum(g1 * log.SE33)
      
      if(sum(g1) != 12 || any(g1 < 0 || g1 > 12))
        browser()
      
      if(l1 > l0 || log.u[i] < l1 - l0){
        g0 = g1
        l0 = l1
        n0 = n0 + 1
      }
    }else if(toss[i] >= 0.5 & n0 > 1){
      #death - pick a
      #cat(paste(i, n0, l0, "death\n"))
      g1 = g0
      SE33subs = SE33
      SE33subs[g1 == 0] = 0
      death = sample(1:52, 1, prob = SE33subs)
      SE33subs[death] = 0
      k = sample(1:52, g1[death], replace = T, prob = SE33subs)
      for(a in k) g1[a] = g1[a] + 1
      g1[death] = 0
      l1 = sum(g1 * log.SE33)
      
      if(sum(g1) != 12 || any(g1 < 0 || g1 > 12))
        browser()
      
      if(l1 > l0 || log.u[i] < l1 - l0){
        g0 = g1
        l0 = l1
        n0 = n0 - 1
      }
    }
    
    res[i] = n0
    setTxtProgressBar(pb, i)
    
  }
  cat(paste(l0, "\n"))
  return(res)
}

