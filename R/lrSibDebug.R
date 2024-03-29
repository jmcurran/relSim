#' Likelihood Ratio / Kinship Index for full-siblings
#' 
#' Calculates Likelihood Ratio comparing the probability of two profiles if
#' they are indeed full-sibs compared to unrelated. This is sometimes called
#' the kinship index (KI) for full-sibs. This function is identical to lrSib
#' except that the calculation is performed in R, and provides full calculation
#' detail at each locus. It exists primarily to check that the correct formula
#' and logic is being applied in the LR calculation so that the result can be
#' manually verified.
#' 
#' 
#' @param sib1 A matrix consisting of 2 columns and nLoci rows. Each entry in
#' the matrix is the (coded) allele held by the individual. This represents the
#' alleged sibling. The relationship is reflexive so it does not matter which
#' profile is labelled sib1 and sib2.
#' @param sib2 See \code{sib1}
#' @param Freqs A list containing two lists labelled loci and freqs. The second
#' list is a list of vectors containing the allele frequencies of each allele
#' at each locus in the multiplex.
#' @return A list containing three elements Lines, lr, and Cases. Lines is a
#' list of strings containing the calculation at each locus so that the result
#' can be written to file for example. Cases is a numeric code listing which
#' logical case (1-11) the locus falls into for the profiles in question. lr is
#' the KI for full-sibs for the two profiles.
#' @author James M. Curran
#' @seealso lrSib, lrPC, IBS
#' @references Buckleton, J, Triggs, C.M., and Walsh, S.J. (2005)\emph{Forensic
#' DNA Evidence Interpretation}, CRC Press., Boca Raton, FL. p.411
#' @examples
#' 
#' data(fbiCaucs)
#' P1 = randomProfile(fbiCaucs)
#' S1 = randomSib(P1, fbiCaucs)
#' P2 = randomProfile(fbiCaucs)
#' cat(paste(lrSibDebug(P1, S1, fbiCaucs)$Lines))
#' cat(paste(lrSibDebug(P1, P2, fbiCaucs)$Lines))
#' 
#' @export lrSibDebug
lrSibDebug = function(sib1, sib2, Freqs){
    lr = 1
    nLoci = length(sib1)/2
    Lines = rep("", nLoci+1)
    Cases = rep(0, nLoci)

    line = sprintf("%8s,%10s,%10s,%6s,%6s,%4s,%28s,%10s,%10s,%6s,%10s\n", "Locus", "Sib 2", "Sib 1", "Rule 1", "Rule 2", "Case","Formula",
                                                                "Pr(A)","Pr(B)", "Value", "Cum. LR")

    Lines[1] = line
    ##cat(line)
    ##line = paste(paste(rep('-',100),collapse=''),"\n")
    ##cat(line)

    for(nLoc in 1:nLoci){
        f = Freqs$freqs[[nLoc]]
        a = names(f)
        i1 = 2*nLoc - 1
        i2 = i1 + 1

        g1 = paste(a[sib1[i1:i2]],collapse="/")
        g2 = paste(a[sib2[i1:i2]],collapse="/")
        r1 = ""
        r2 = ""
        pA = 0
        pB = 0
        m = ""
        mn = 0
        case = 0

        if(sib2[i1] == sib2[i2]){ ## sib2 is aa
            pA = f[sib2[i1]]
            r2 = "aa"

            if(sib1[i1] == sib1[i2]){ ## sib1 is aa or bb
                if(sib1[i1] == sib2[i1]){ ## sib1 is aa
                    r1 = "aa"
                    m = "(1+pA)*(1+pA)/(4*pA*pA)"
                    mn = (1+pA)*(1+pA)/(4*pA*pA)
                    case = 1
                    lr = lr * mn
                }else{ ## sib1 is bb
                    r1 = "bb"
                    m = "1/4"
                    mn = 0.25
                    case = 2
                    lr = lr * mn
                }
            }else{ ## sib2 is ab or bc
                if(sib1[i1] != sib2[i1] & sib1[i2] != sib2[i1]){ ## sib1 is bc
                    r1 = "bc"
                    m = "1/4"
                    mn = 0.25
                    case = 3
                    lr = lr * mn
                }else{ ## sib1 is ab
                    r1 = "ab"
                    m = "(1+pA)/(4*pA)"
                    mn = (1+pA)/(4*pA)
                    case = 4
                    lr = lr * mn
                }
            }
        }else{ ## sib2 is ab
            pA = f[sib2[i1]]
            pB = f[sib2[i2]]

            r2 = "ab"

            if(sib1[i1] == sib1[i2]){ ## sib1 is aa, bb or cc
                if(sib1[i1] == sib2[i1]){ ## sib1 is aa
                    r1 = "aa"
                    m = "(1+pA)/(4*pA)"
                    mn = (1+pA)/(4*pA)
                    case = 5
                    lr = lr * mn
                }else if(sib1[i1] == sib2[i2]){ ### sib1 is bb
                    r1 = "bb"
                    m = "(1+pB)/(4*pB)"
                    mn = (1+pB)/(4*pB)
                    case = 6
                    lr = lr * mn
                }else{ ## sib1 is cc
                    r1 = "cc"
                    m = "1/4"
                    mn = 0.25
                    case = 7
                    lr = lr * mn
                }
            }else{ ## sib1 is ab, ac, bc, or cd
                if(sib1[i1] != sib2[i1] & sib1[i2] != sib2[i2]
                   & sib1[i2] != sib2[i1] & sib1[i1] != sib2[i2]){ ## sib1 is cd
                    r1 = "cd"
                    m = "1/4"
                    mn = 0.25
                    case = 8
                    lr = lr * mn
                }else if(sib1[i1] == sib2[i1] & sib1[i2] == sib2[i2]){ ## sib1 is ab
                    r1 = "ab"
                    m = "(1+pA+pB+2*pA*pB)/(8*pA*pB)"
                    mn = (1+pA+pB+2*pA*pB)/(8*pA*pB)
                    case = 9
                    lr = lr * mn
                }else{ ## sib1 is ac or bc
                    if((sib1[i1] == sib2[i1] & sib1[i2] != sib2[i2])
                       | (sib1[i2] == sib2[i1] & sib1[i1] !=sib2[i2])){ ## sib1 is ac
                        r1 = "ac"
                        m = "(1+2*pA)/(8*pA)"
                        mn = (1+2*pA)/(8*pA)
                        case =10
                        lr = lr * mn
                    }else{ ## sib1 is bc
                        r1 = "bc"
                        m = "(1+2*pB)/(8*pB)"
                        mn = (1+2*pB)/(8*pB)
                        case = 11
                        lr = lr * mn
                    }
                }
            }
        }

        loc = Freqs$loci[nLoc]
        line = sprintf("%8s,%10s,%10s,%6s,%6s,%4d,%28s,%10.8f,%10.8f,%6.4f,%10.2f\n", loc, g2, g1, r2, r1, case, m, pA, pB, mn, lr)
        ##at(line)
        Lines[nLoc + 1] = line
        Cases[nLoc] = case
    }

    invisible(list(Lines = Lines, lr = lr, Cases = Cases))
}
