randomProfiles = function(Freqs, BlockSize = 1000){
    Profile = vector(mode = "list", length = BlockSize)
    profileVec = .randomProfiles(Freqs$freq, BlockSize)
    nLoci = nLoci = length(Freqs$loci)
    
    for(b in 1:BlockSize){
        i1 = (b - 1)*2*nLoci + 1
        i2 =  b*2*nLoci

        Profile[[b]] = matrix(profileVec[i1:i2], ncol = 2, nrow = nLoci, byrow = T)
        class(Profile[[b]]) = "profile"
    }

    return(Profile)
}

# randomDupProfiles = function(Freqs, N = 1000, numFounders = 5){
#   Profile = vector(mode = "list", length = N)
#   profileVec = .randomProfiles(Freqs$freq, numFounders)
#   nLoci = nLoci = length(Freqs$loci)
#   
#   i = sort(sample(1:numFounders, size = N, replace = TRUE))
#   
#   for(j in seq_along(i)){
#     i1 = (i[j] - 1) * 2 * nLoci + 1
#     i2 =  i[j] * 2 * nLoci
#     
#     Profile[[j]] = matrix(profileVec[i1:i2], ncol = 2, nrow = nLoci, byrow = T)
#     class(Profile[[j]]) = "profile"
#   }
#   
#   return(Profile)
# }
# 
# mutateProfiles = function(Profiles, freqs, maxMutatations, numMutations, mutationRate = 0){
#   
#   N = length(Profiles)
#   
#   if(mutationRate != 0){
#     if(mutationRate < 0 || mutationRate > 1){
#       stop("Mutation rate should be between 0 and 1")
#     }
#     numMutations = N * mutationRate
#   }
#   
#   i = sort(sample.int(N, numMutations))
#   numMutations = sample.int(maxMutatations, numMutations, replace = TRUE)
#   
#   for(j in i){
#     p = Profiles[[j]]
#     pos = sample.int(1:(2 * nLoci))
#   }
# }
