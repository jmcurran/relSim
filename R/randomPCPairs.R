randomPCPairs = function(Freqs, BlockSize = 1){
    ## nLoci = length(Freqs$loci)
    ## profChild = matrix(0, nc = 2, nr = nLoci)

    ## for(nLoc in 1:nLoci){
    ##     f = Freqs$freqs[[nLoc]]
    ##     a = sample(1:length(f), 1, prob = f)
    ##     u = runif(1)

    ##     if(u < 0.5){
    ##         profChild[nLoc,] = c(profile[nLoc,1], a)
    ##     }else{
    ##         profChild[nLoc,] = c(a, profile[nLoc,2])
    ##     }

    ##     if(profChild[nLoc,1] > profChild[nLoc, 2]){
    ##         swap = profChild[nLoc, 1]
    ##         profChild[nLoc, 1] = profChild[nLoc, 2]
    ##         profChild[nLoc, 2] = swap
    ##     }
    ## }

    ## return(profChild)

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
    pVecParent = rep(0, 2*nLoci*BlockSize)
    pVecChild = rep(0, 2*nLoci*BlockSize)


    pVecParent = .C("randomProfiles", p = as.integer(pVecParent),
                               nLoci = as.integer(nLoci),
                               f = as.double(f), n = as.integer(n),
                               u = as.double(runif(2*nLoci*BlockSize)),
                               b = as.integer(BlockSize))$p

    pVecChild = .C("randomChildren", pParent = as.integer(pVecParent),
                               pChild = as.integer(pVecChild),
                               nLoci = as.integer(nLoci),
                               f = as.double(f), n = as.integer(n),
                               u = as.double(runif(3*nLoci*BlockSize)),
                               b = as.integer(BlockSize))$pChild



    for(b in 1:BlockSize){
        i1 = (b - 1)*2*nLoci + 1
        i2 =  b*2*nLoci

        Profile[[b]] = vector(mode = "list", length = 2)
        names(Profile[[b]]) = c("parent", "child")

        Profile[[b]]$parent = pVecParent[i1:i2]
        class(Profile[[b]]$parent) = "profile"

        Profile[[b]]$child = pVecChild[i1:i2]
        class(Profile[[b]]$child) = "profile"
    }

    if(BlockSize==1){
        return(Profile[[1]])
    }else{
        return(Profile)
    }
}
