lrSib = function(sib1, sib2, Freqs  = NULL, nLoci = length(sib1)/2,
                 f = NULL, n = NULL){

    if(!is.null(Freqs)){
        f = unlist(Freqs$freqs)
        n = sapply(Freqs$freqs, length)
    }

    if(is.null(f) | is.null(n)){
        stop("Must specify either Freqs or f and n")
    }


    lr = 1
    lr = .C("lrSib", sib1 = as.integer(sib1), sib2 = as.integer(sib2),
                     nLoci = as.integer(nLoci), f = as.double(f),
                     n = as.integer(n), lr = as.double(lr))$lr
    ## for(nLoc in 1:nLoci){
    ##     f = Freqs$freqs[[nLoc]]
    ##     llr = 0
    ##     i1 = 2*nLoc - 1
    ##     i2 = i1 + 1

    ##     lr = lr * .C("locusLRSib", sib1 = as.integer(c(sib1[i1:i2])),
    ##                                sib2 = as.integer(c(sib2[i1:i2])),
    ##                                f = as.double(f), llr = as.double(llr))$llr

        ## if(sib2[nLoc,1] == sib2[nLoc,2]){ ## sib2 is aa
    ##         pA = f[sib2[nLoc,1]]

    ##         if(sib1[nLoc,1] == sib1[nLoc,2]){ ## sib1 is aa or bb
    ##             if(sib1[nLoc,1] == sib2[nLoc, 1]){ ## sib1 is aa
    ##                 lr = lr * (1+pA)*(1+pA)/(4*pA*pA)
    ##             }else{ ## sib1 is bb
    ##                 lr = lr * 0.25
    ##             }
    ##         }else{ ## sib2 is ab or bc
    ##             if(sib1[nLoc,1] != sib2[nLoc,1] & sib1[nLoc, 2] != sib2[nLoc, 1]){ ## sib1 is bc
    ##                 lr = lr * 0.25
    ##             }else{ ## sib1 is ab
    ##                 lr = lr * (1+pA)/(4*pA)
    ##             }
    ##         }
    ##     }else{ ## sib2 is ab
    ##         pA = f[sib2[nLoc,1]]
    ##         pB = f[sib2[nLoc,2]]

    ##         if(sib1[nLoc, 1] == sib1[nLoc, 2]){ ## sib1 is aa, bb or cc
    ##             if(sib1[nLoc,1] == sib2[nLoc,1]){ ## sib1 is aa
    ##                 lr = lr * (1+pA)/(4*pA)
    ##             }else if(sib1[nLoc,1] == sib2[nLoc,2]){ ### sib1 is bb
    ##                 lr = lr * (1+pB)/(4*pB)
    ##             }else{ ## sib1 is cc
    ##                 lr = lr * 0.25
    ##             }
    ##         }else{ ## sib1 is ab, ac, bc, or cd
    ##             if(sib1[nLoc,1] != sib2[nLoc,1] & sib1[nLoc,2] != sib2[nLoc,2]
    ##                & sib1[nLoc,2] != sib2[nLoc,1] & sib1[nLoc,1] != sib2[nLoc,2]){ ## sib1 is cd
    ##                 lr = lr * 0.25
    ##             }else if(sib1[nLoc,1] == sib2[nLoc,1] & sib1[nLoc,2] == sib2[nLoc,2]){ ## sib1 is ab
    ##                 lr = lr * (1+pA+pB+2*pA*pB)/(8*pA*pB)
    ##             }else{ ## sib1 is ac or bc
    ##                 if((sib1[nLoc, 1] == sib2[nLoc, 1] & sib1[nLoc,2] != sib2[nLoc,2])
    ##                    | (sib1[nLoc, 2] == sib2[nLoc,1] & sib1[nLoc,1] !=sib2[nLoc, 2])){ ## sib1 is ac
    ##                     lr = lr * (1+2*pA)/(8*pA)
    ##                 }else{ ## sib1 is bc
    ##                     lr = lr * (1+2*pB)/(8*pB)
    ##                 }
    ##             }
    ##         }
    ##     }
    ##}

    return (lr)
}
