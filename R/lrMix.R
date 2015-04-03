lrMix = function(profiles, Freqs){
    N  = length(profiles)
    results = matrix(0, nrow = N, ncol = nLoci)

    for(i in 1:N){
      results[i,] = LRmix_Caller(profiles[[i]][[1]],
                        profiles[[i]][[2]],
                        Freqs$freqs)
    }

    colnames(results) = Freqs$loci
    return(results)
}
