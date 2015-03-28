readFreqs = function(strPath){
    f1 = file(strPath, "r")

    if(isOpen(f1)){
        Lines = readLines(f1)
        nLines = length(Lines)
        ## cat(paste("Read", nLines, "lines from", strPath,"\n"))
        close(f1)
    }else{
        stop(paste("Couldn't open :", strPath, "for reading"))
    }

    Freqs = NULL

    nLoci = as.numeric(Lines[1])
    ## cat(paste(nLoci, "\n"))

    Loci = rep("", nLoci)
    freqs = vector(nLoci, mode = "list")
    currLine = 2

    for(nLoc in 1:nLoci){
        Tokens = unlist(strsplit(Lines[currLine], ","))
        currLine = currLine + 1

        Loci[nLoc] = Tokens[1]
        nAlleles = as.numeric(Tokens[2])
        Alleles = rep("", nAlleles)

        freqs[[nLoc]] = rep(0, nAlleles)

        for(nA in 1:nAlleles){
            Tokens = unlist(strsplit(Lines[currLine], ","))
            currLine = currLine + 1
            Alleles[nA] = Tokens[2]
            freqs[[nLoc]][nA] = as.numeric(Tokens[3])
        }

        names(freqs[[nLoc]]) = Alleles
    }

    names(freqs) = Loci
    Freqs = list(loci = Loci, freqs = freqs)

    return(Freqs)
}
