writePop = function(fileName,  pop, addAmelo = FALSE, delim = ','){
  nLoci = pop$nLoci
  Alleles = lapply(pop$Freqs$freqs, names)
  Alleles = lapply(Alleles, function(x){x[grep("R", x)] = "108.1"; x})
  
  toProf = function(prof){
    A = rep("", 2 * nLoci)
    if(addAmelo){
      A = c(A, "1", "2")
    }
    
    for(i in 1:nLoci){
      i1 = 2 * i - 1
      i2 = i1 + 1
      A[i1] = Alleles[[i]][prof[i1]]
      A[i2] = Alleles[[i]][prof[i2]]
    }
    
    return(paste(A, collapse = delim))
  }
  
  popMatrix = matrix(pop$profiles, nrow = pop$nProfiles, byrow = T) 
  profiles = apply(popMatrix, 1, toProf)
  f1 = file(fileName, 'w')
  locusHeader = rep(pop$Freqs$loci, rep(2, nLoci))
  if(addAmelo)
    locusHeader = paste0(locusHeader, ',"Amelo", "Amelo"')
  writeLines(paste0('"', paste(rep(pop$Freqs$loci, rep(2, nLoci)), collapse = '","'), '"'), f1)
  writeLines(profiles, f1)
  close(f1)
}
