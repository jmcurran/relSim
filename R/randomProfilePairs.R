randomProfilePairs = function(Freqs, BlockSize = 1){
    ## nLoci = length(Freqs$loci)
    ## profile = matrix(0, nc = 2, nr = nLoci)

    ## for(nLoc in 1:nLoci){
    ##     f = Freqs$freqs[[nLoc]]
    ##     profile[nLoc,] = sample(1:length(f), 2, replace = T, prob = f)

    ##     if(profile[nLoc,1] > profile[nLoc, 2]){
    ##         swap = profile[nLoc, 1]
    ##         profile[nLoc, 1] = profile[nLoc, 2]
    ##         profile[nLoc, 2] = swap
    ##     }
    ## }

    ## class(profile) = "profile"
    ## return(profile)

    f = unlist(Freqs$freqs)
    n = sapply(Freqs$freqs, length)
    nLoci = length(Freqs$loci)


    Profile = vector(mode = "list", length = BlockSize)
    pVec1 = rep(0, 2*nLoci*BlockSize)
    pVec2 = rep(0, 2*nLoci*BlockSize)

    u1 = runif(2*nLoci*BlockSize)
    pVec1 = .C("randomProfiles", p = as.integer(pVec1),
                               nLoci = as.integer(nLoci),
                               f = as.double(f), n = as.integer(n),
                               u = as.double(u1),
                               b = as.integer(BlockSize))$p

    u1 = runif(2*nLoci*BlockSize)
    pVec2 = .C("randomProfiles", p = as.integer(pVec2),
                               nLoci = as.integer(nLoci),
                               f = as.double(f), n = as.integer(n),
                               u = as.double(u1),
                               b = as.integer(BlockSize))$p

    for(b in 1:BlockSize){
        i1 = (b - 1)*2*nLoci + 1
        i2 =  b*2*nLoci

        Profile[[b]]$prof1 = pVec1[i1:i2]
        Profile[[b]]$prof2 = pVec2[i1:i2]
        class(Profile[[b]]$prof1) = "profile"
        class(Profile[[b]]$prof2) = "profile"
    }

    if(BlockSize==1){
        return(Profile[[1]])
    }else{
        return(Profile)
    }
}
