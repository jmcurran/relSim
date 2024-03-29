#' Export a population with substructure to a Nexus file
#' 
#' Exports a population with population substructure to a Nexus formatted file
#' so that GDA can be used to check the Fst calculations
#' 
#' 
#' @param Pop An object of type 'population' - see \code{breedFst} for a
#' description of the object
#' @param fileName The name of the file output file
#' @author James M. Curran
#' @seealso breedFst
#' @references Maddison DR, Swofford DL, Maddison WP (1997), NEXUS: An
#' extensible file format for systematic information, Systematic Biology 46
#' (4): 590--621.
#' 
#' Zaykin, D. and Lewis, P., GDA - software to accompany Genetic Data Analysis
#' II, <http://phylogeny.uconn.edu/software/>.
#' @examples
#' 
#' ## Don't run
#' \dontrun{
#' data(USCaucs)
#' p = breedFst(USCaucs)
#' toNexus(p)
#' }
#' 
#' @export toNexus
toNexus = function(Pop, fileName = 'output.nex'){
    if(!is(Pop, "population"))
        stop("Pop must be an object of class 'population'")

    nLoci = Pop$nLoci
    ns = Pop$nSubpops

    f1 = file(fileName, 'w')

    if(!isOpen(f1)){
        stop(paste("Could not open", fileName, "for writing\n"))
    }

    writeLines(c("#nexus", "", "begin gdadata;"), f1)

    line = paste("\tdimensions nloci=", nLoci, "npops=", ns, ";")
    writeLines(line, f1)

    writeLines(c("\tformat tokens labels missing=? datapoint=standard;",
                 "\tlocusallelelabels"), f1);

    locusLines = paste('\t\t', 1:nLoci, ' ',
                       "'", paste('Locus_', Pop$Freqs$loci, sep = ''),
                       "',", sep = "")
    locusLines[nLoci] = substr(locusLines[nLoci], 0,
                               nchar(locusLines[nLoci]) - 1)

    writeLines(locusLines, f1)
    writeLines(c(';', "MATRIX"), f1)

    pop = matrix(Pop$profiles, ncol = 2 * nLoci, byrow = T)
    Ns = nrow(pop) / ns

    toProf = function(x){
        n = length(x) / 2;
        i1 = 2 * (1:n) - 1; i2 = i1 + 1;
        paste(paste(x[i1], '/', x[i2], sep = ""), collapse = " ")
    }

    for(s in 1:ns){
        popLabel = paste('Pop_', s, ':', sep = '')
        writeLines(popLabel, f1)

        i1 = (s - 1)*Ns + 1
        i2 = s*Ns

        writeLines(paste('\t\t', 1:Ns, apply(pop[i1:i2,], 1, toProf)), f1)

        if(s <= (ns - 1))
            writeLines(',', f1)
    }
    writeLines(c("\t;", "END;"), f1)
    close(f1)

}
