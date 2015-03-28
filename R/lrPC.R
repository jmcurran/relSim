lrPC = function(parent, child, Freqs = NULL,
                nLoci = length(parent)/2, f = NULL, n = NULL){

    if(!is.null(Freqs)){
        f = unlist(Freqs$freqs)
        n = sapply(Freqs$freqs, length)
    }

    if(is.null(f) | is.null(n)){
        stop("Must specify either Freqs or f and n")
    }

    lr = 1
    lr = .C("lrPC",  parent = as.integer(parent), child = as.integer(child),
                     nLoci = as.integer(nLoci), f = as.double(f),
                     n = as.integer(n), lr = as.double(lr))$lr


    ## lr = 1
    ## nLoci = length(parent)/2

    ## for(nLoc in 1:nLoci){
    ##     f = Freqs$freqs[[nLoc]]
    ##     llr = 1
    ##     i1 = 2*nLoc - 1
    ##     i2 = i1 + 1

    ##     llr = .C("locusLRPC", parent = as.integer(parent[i1:i2]),
    ##                           child = as.integer(child[i1:i2]),
    ##                           f = as.double(f),
    ##                           llr = as.double(llr))$llr

    ##     if(llr == 0)
    ##         return(0)
    ##     else
    ##         lr = lr*llr

        ## if(child[nLoc,1] == child[nLoc,2]){ ## child is aa
        ##     pA = f[child[nLoc,1]]

        ##     if(parent[nLoc,1] == parent[nLoc,2]){ ## parent is aa or bb
        ##         if(parent[nLoc,1] == child[nLoc,1]){ ## parent is aa
        ##             lr = lr / pA
        ##         }else{ ## parent is bb
        ##             return(0)
        ##         }
        ##     }else{ ## parent is ab or bc
        ##         if(parent[nLoc,1] != child[nLoc,1] & parent[nLoc,2] != child[nLoc, 1]){ ## parent is bc
        ##             return(0)
        ##         }else{ ## parent is ab
        ##             lr = lr / (2*pA)
        ##         }
        ##     }
        ## }else{ ## child is ab
        ##     pA = f[child[nLoc,1]]
        ##     pB = f[child[nLoc,2]]

        ##     if(parent[nLoc,1] == parent[nLoc,2]){ ## parent is aa, bb, or cc
        ##         if(parent[nLoc,1] == child[nLoc,1]){ ## parent is aa
        ##             lr = lr / (2*pA)
        ##         }else if(parent[nLoc,1] == child[nLoc,2]){ ## parent is bb
        ##             lr = lr / (2*pB)
        ##         }else{ ## parent is cc
        ##             return(0);
        ##         }
        ##     }else{ ## parent is ab, ac, bc, or cd
        ##         if(parent[nLoc,1] == child[nLoc, 1] & parent[nLoc,2] == child[nLoc, 2]){ ## parent is ab
        ##             lr = lr * (pA+pB)/(4*pA*pB)
        ##         }else if(parent[nLoc,1] == child[nLoc,1] | parent[nLoc,2] == child[nLoc, 1]){ ## parent is ac or ca
        ##             lr = lr / (4*pA)
        ##         }else if(parent[nLoc,1] == child[nLoc,2] | parent[nLoc,2] == child[nLoc, 2]){ ## parent is bc or cb
        ##             lr = lr /(4*pB)
        ##         }else{ ## parent is cd
        ##             return(0)
        ##         }
        ##     }
        ## }
    ##}
    return(lr)
}
