calcFst = function(Pop, subPopIdx = NULL){
    if(class(Pop) != "population")
        stop("Pop must be of class Population")

    nLoci = Pop$nLoci
    NumLocusAlleles = sapply(Pop$Freqs$freqs, length)

    if(is.data.frame(Pop$profiles)){
        Pop$profiles = as.vector(t(as.matrix(Pop$profiles[,-1])))
    }

    N = length(Pop$profiles)

    if(N %% nLoci !=0){
        stop("This proceedure only works for complete profiles")
    }

    N = Pop$nProfiles
    ns = Pop$nSubpops

    if(N %% ns != 0)
        stop("The number of subpopulations must evenly divide N")

    if(is.null(subPopIdx)){
        Ns = N / ns
        subPopIdx = rep(1:ns, rep(Ns, ns))
    }else{
        if(length(subPopIdx != N))
            stop("subPopIdx length must equal pop length")

        ns = length(unique(subPopIdx))
    }

    Fst = rep(0, nLoci + 1)

    r = .C("calcFst", profiles = as.integer(Pop$profiles),
                      subPopIdx = as.integer(subPopIdx),
                      N = as.integer(N), ns = as.integer(ns),
                      nLoci = as.integer(nLoci),
                      NumLocusAlleles = as.integer(NumLocusAlleles),
                      Fst = as.double(Fst))


    Fst = r$Fst
    names(Fst)[1:nLoci] = names(Pop$loci)
    names(Fst)[nLoci + 1] = "Overall"
    return(Fst)
}
