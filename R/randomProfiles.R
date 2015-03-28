randomProfiles = function(Freqs, BlockSize = 1000){
    f = unlist(Freqs$freqs)
    n = sapply(Freqs$freqs, length)
    nLoci = length(Freqs$loci)


    Profile = vector(mode = "list", length = BlockSize)
    pVec = rep(0, 2*nLoci*BlockSize)

    pVec = .C("randomProfiles", p = as.integer(pVec), nLoci = as.integer(nLoci),
                               f = as.double(f), n = as.integer(n),
                               u = as.double(runif(2*nLoci*BlockSize)),
                               b = as.integer(BlockSize))$p

    for(b in 1:BlockSize){
        i1 = (b - 1)*2*nLoci + 1
        i2 =  b*2*nLoci

        Profile[[b]] = matrix(pVec[i1:i2], ncol = 2, nrow = nLoci, byrow = T)
        class(Profile[[b]]) = "profile"
    }

    return(Profile)
}
