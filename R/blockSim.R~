blockSim = function(N, Freqs, rel = "UN", ibsthresh = NULL, kithresh = NULL,
                    code = 1, falseNeg = TRUE, BlockSize = N/10){

    rel = toupper(rel)
    if(!grepl("(UN|FS|PC)", rel)){
        stop("rel must be one of 'UN', 'FS' or 'PC'")
    }

    f = unlist(Freqs$freqs)
    n = sapply(Freqs$freqs, length)
    nLoci = length(Freqs$loci)


    pVec1 = rep(0, 2*nLoci*BlockSize)
    pVec2 = rep(0, 2*nLoci*BlockSize)

    nBlocks = N/BlockSize
    if(is.null(ibsthresh) & is.null(kithresh))
        stop("You must specify one or both of ibsthresh or kithresh")

    nResults = 0
    if(is.null(ibsthresh)){
        nResults = length(kithresh)
        ibsthresh = rep(0, nResults) ## dummy vals
    }else{
        nResults = length(ibsthresh)
        kithresh = rep(0, nResults)
    }

    if(nResults == 0)
        stop("Noththing to count")

    nTotal = rep(0, nResults)

    pb = txtProgressBar(min = 0, max = nBlocks, style = 3)

    if(rel == "UN"){
        for(block in 1:nBlocks){
            pVec1 = .C("randomProfiles", p = as.integer(pVec1),
                        nLoci = as.integer(nLoci),
                        f = as.double(f), n = as.integer(n),
                        u = as.double(runif(2*nLoci*BlockSize)),
                        b = as.integer(BlockSize))$p

            pVec2 = .C("randomProfiles", p = as.integer(pVec2),
                        nLoci = as.integer(nLoci),
                        f = as.double(f), n = as.integer(n),
                        u = as.double(runif(2*nLoci*BlockSize)),
                        b = as.integer(BlockSize))$p

            count = rep(0, length(ibsthresh))
            nc = length(count)
            count = .C("blockStatCounts", as.integer(pVec1), as.integer(pVec2),
                        as.integer(nLoci), as.integer(BlockSize),
                        f = as.double(f), n = as.integer(n),
                        code = as.integer(code),
                        falseNeg = as.integer(falseNeg),
                        ibs = as.integer(ibsthresh),
                        ki = as.double(kithresh),
                        count = as.integer(count),
                        nc = as.integer(nc))$count

            nTotal = nTotal + count
            setTxtProgressBar(pb, block)
        }
    }else if(rel == "FS"){
        for(block in 1:nBlocks){
            pVec1 = .C("randomProfiles", prof1 = as.integer(pVec1),
                        nLoci = as.integer(nLoci),
                        f = as.double(f), n = as.integer(n),
                        u = as.double(runif(2*nLoci*BlockSize)),
                        b = as.integer(BlockSize))$prof1

            pVec2 = .C("randomSibs", prof1 = as.integer(pVec1),
                        prof2 = as.integer(pVec2),
                        nLoci = as.integer(nLoci),
                        f = as.double(f), n = as.integer(n),
                        u = as.double(runif(3*nLoci*BlockSize)),
                        b = as.integer(BlockSize))$prof2

            count = rep(0, length(ibsthresh))
            nc = length(count)
            count = .C("blockStatCounts", as.integer(pVec1), as.integer(pVec2),
                        as.integer(nLoci), as.integer(BlockSize),
                        f = as.double(f), n = as.integer(n),
                        code = as.integer(code),
                        falseNeg = as.integer(falseNeg),
                        ibs = as.integer(ibsthresh),
                        ki = as.double(kithresh),
                        count = as.integer(count),
                        nc = as.integer(nc))$count

            nTotal = nTotal + count
            setTxtProgressBar(pb, block)
        }
    }else if(rel == "PC"){
        for(block in 1:nBlocks){
            pVec1 = .C("randomProfiles", prof1 = as.integer(pVec1),
                        nLoci = as.integer(nLoci),
                        f = as.double(f), n = as.integer(n),
                        u = as.double(runif(2*nLoci*BlockSize)),
                        b = as.integer(BlockSize))$prof1

            pVec2 = .C("randomChildren", prof1 = as.integer(pVec1),
                        prof2 = as.integer(pVec2),
                        nLoci = as.integer(nLoci),
                        f = as.double(f), n = as.integer(n),
                        u = as.double(runif(2*nLoci*BlockSize)),
                        b = as.integer(BlockSize))$prof2

            count = rep(0, length(ibsthresh))
            nc = length(count)
            count = .C("blockStatCounts", as.integer(pVec1), as.integer(pVec2),
                        as.integer(nLoci), as.integer(BlockSize),
                        f = as.double(f), n = as.integer(n),
                        code = as.integer(code),
                        falseNeg = as.integer(falseNeg),
                        ibs = as.integer(ibsthresh),
                        ki = as.double(kithresh),
                        count = as.integer(count),
                        nc = as.integer(nc))$count

            nTotal = nTotal + count
            setTxtProgressBar(pb, block)
        }
    }

    return(list(nTotal = nTotal, N = N, p = nTotal/N))
}
