lrMix = function(profiles, Freqs){
    n = c(0,cumsum(sapply(Freqs$freqs, length)))
    f = unlist(Freqs$freqs)
    nLoci = length(Freqs$loci)
    N  = length(profiles)
    results = matrix(0, nrow = N, ncol = nLoci)

    for(i in 1:N){
        r = .C("LRmix", v = as.integer(profiles[[i]][[1]]),
                        s = as.integer(profiles[[i]][[2]]),
                        nLoc = as.integer(nLoci), f = as.double(f),
                        offset = as.integer(n),
                        result = as.double(results[i,]))
        results[i,] = r$result
    }

    colnames(results) = Freqs$loci
    return(results)
}
