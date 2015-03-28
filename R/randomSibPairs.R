randomSibPairs = function(Freqs, BlockSize = 1){
    ## nLoci = length(Freqs$loci)
    ## profSib = matrix(0, nc = 2, nr = nLoci)

    ## for(nLoc in 1:nLoci){
    ##     f = Freqs$freqs[[nLoc]]
    ##     i = sample(1:4, 1)
    ##     a = sample(1:length(f), 2, replace = TRUE, prob = f)

    ##     switch(i,
    ##            {profSib[nLoc,] = profile[nLoc,]},
    ##            {profSib[nLoc,] = c(profile[nLoc,1], a[1])},
    ##            {profSib[nLoc,] = c(a[1], profile[nLoc,2])},
    ##            {profSib[nLoc,] = a}
    ##            )


    ##     if(profSib[nLoc,1] > profSib[nLoc, 2]){
    ##         swap = profSib[nLoc, 1]
    ##         profSib[nLoc, 1] = profSib[nLoc, 2]
    ##         profSib[nLoc, 2] = swap
    ##     }
    ## }

    f = unlist(Freqs$freqs)
    n = sapply(Freqs$freqs, length)
    nLoci = length(Freqs$loci)

    Profile = vector(mode = "list", length = BlockSize)
    pVecSib1 = rep(0, 2*nLoci*BlockSize)
    pVecSib2 = rep(0, 2*nLoci*BlockSize)

    pVecSib1 = .C("randomProfiles", p = as.integer(pVecSib1),
                               nLoci = as.integer(nLoci),
                               f = as.double(f), n = as.integer(n),
                               u = as.double(runif(2*nLoci*BlockSize)),
                               b = as.integer(BlockSize))$p

    pVecSib2 = .C("randomSibs", pSib1 = as.integer(pVecSib1),
                               pSib2 = as.integer(pVecSib2),
                               nLoci = as.integer(nLoci),
                               f = as.double(f), n = as.integer(n),
                               u = as.double(runif(3*nLoci*BlockSize)),
                               b = as.integer(BlockSize))$pSib2



    for(b in 1:BlockSize){
        i1 = (b - 1)*2*nLoci + 1
        i2 =  b*2*nLoci

        Profile[[b]] = vector(mode = "list", length = 2)
        names(Profile[[b]]) = c("sib1", "sib2")

        Profile[[b]]$sib1 = pVecSib1[i1:i2]
        class(Profile[[b]]$sib1) = "profile"

        Profile[[b]]$sib2 = pVecSib2[i1:i2]
        class(Profile[[b]]$sib2) = "profile"
    }

    if(BlockSize==1){
        return(Profile[[1]])
    }else{
        return(Profile)
    }
}
