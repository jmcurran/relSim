breedFst = function(Freqs, theta = 0.01, N = 10000, ns = 10,
                    DNAtools = FALSE){
    if(N<1000){
        stop("N must be >= 1000")
    }

    if(N %% ns != 0){
        stop("ns must divide N into a whole number\tThat is the subpopulation sizes must be integers")
    }

    if(theta <=0 || theta >= 0.5){
        stop("0 < theta < 0.5")
    }

    Ns = N / ns

    nGen = ceiling(log(1 - theta) / log(1 - 1/(2*Ns)))
    cat(paste("Breeding for", nGen, "generations\n"))

    f = unlist(Freqs$freqs)
    n = sapply(Freqs$freqs, length)
    nLoci = length(Freqs$freqs)

    ## generate the parental population

    pVec = rep(0, 2 * nLoci * N)
    pVec = .C("randomProfiles", p = as.integer(pVec), nLoci = as.integer(nLoci),
                               f = as.double(f), n = as.integer(n),
                               u = as.double(runif(2*nLoci*N)),
                               b = as.integer(N))$p


    seeds = floor(29998*runif(3)) + 1
    pb = txtProgressBar(min = 0, max = nGen, style = 3)

    for(t in 1:nGen){
        r = .C("breed", p = as.integer(pVec), ns = as.integer(ns),
                           Ns = as.integer(Ns), nLoci = as.integer(nLoci),
                           x = as.integer(seeds[1]), y = as.integer(seeds[2]),
                           z = as.integer(seeds[3]))
        pVec = r$p
        seeds = with(r, c(x,y,z))
        setTxtProgressBar(pb, t)
    }
    setTxtProgressBar(pb, nGen)

    if(DNAtools){
        pVec = data.frame(cbind(1:N, matrix(pVec, nrow = N, byrow = T)))
        colNames = paste(rep(Freqs$loci, rep(2, nLoci)), c('a','b'), sep = '_')
        names(pVec) = c('id', colNames)
    }

    pop = list(profiles = pVec, nProfiles = N,  nSubpops = ns,
               nLoci = nLoci, theta = theta, Freqs = Freqs)
    class(pop) = "population"

    return(pop)
}
