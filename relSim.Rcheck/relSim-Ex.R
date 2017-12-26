pkgname <- "relSim"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "relSim-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('relSim')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("IBS")
### * IBS

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: IBS
### Title: Identity by state
### Aliases: IBS

### ** Examples


data(fbiCaucs)
P1 = randomProfile(fbiCaucs)
C1 = randomChild(P1, fbiCaucs)
IBS(P1, C1)
IBS(P1, C1, bPrint = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("IBS", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("USCaucs")
### * USCaucs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: USCaucs
### Title: CODIS STR Loci allele frequency data
### Aliases: USCaucs
### Keywords: datasets

### ** Examples


data(USCaucs)
names(USCaucs)
USCaucs$loci
names(USCaucs$freqs)
USCaucs$freqs[[1]]
names(USCaucs$freqs[[1]])
USCaucs$freqs[[1]][1]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("USCaucs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("blockSim")
### * blockSim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: blockSim
### Title: Perform relatives simulations using large memory blocks in C
### Aliases: blockSim

### ** Examples


## not run
## this counts the number of unrelated pairs that are falsely identified
## as siblings using the policy that there are 16 or more matching
## alleles, and the LR/KI is greater than 100,000
## this is a very rare event for the FBI Caucasians with a frequency of
## about 4-5 times in 10 million pairs
## Not run: 
##D data(fbiCaucs)
##D N = 1e8
##D ki = 1e5
##D ibs = 16
##D code = 5
##D BlockSize = 1e6
##D blockSim(N, fbiCaucs, rel = "UN", ibsthresh = ibs, kithresh = ki,
##D          code = code, falseNeg = FALSE, BlockSize = BlockSize)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("blockSim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("breedFst")
### * breedFst

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: breedFst
### Title: Breed a population with an approximate level of theta(Fst)
### Aliases: breedFst

### ** Examples


data(USCaucs)
pop = breedFst(USCaucs)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("breedFst", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcFst")
### * calcFst

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcFst
### Title: Caculate locus-wise and population Fst values
### Aliases: calcFst

### ** Examples


data(USCaucs)
p = breedFst(USCaucs)
fst = calcFst(p)
fst




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcFst", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("checkFreqs")
### * checkFreqs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: checkFreqs
### Title: Make sure that the frequencies are such
### Aliases: checkFreqs

### ** Examples


data(fbiCaucs)
checkFreqs(fbiCaucs)

## induce an error
fbiCaucs$freqs[[1]] = runif(10)
checkFreqs(fbiCaucs)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("checkFreqs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("errorRate")
### * errorRate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: errorRate
### Title: Returns the false positive or false negative rates for a set of
###   IBS and/or KI thresholds
### Aliases: errorRate

### ** Examples


## not run
## Not run: 
##D data(fbiCaucs)
##D unrel = sim(10000)
##D errorRate(unrel)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("errorRate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("exclusionPower")
### * exclusionPower

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: exclusionPower
### Title: Calculate the exclusion power of a multiplex by locus
### Aliases: exclusionPower ep

### ** Examples


data(USCaucs)
ep(USCaucs)

## get the multiplex wide exclusion power
1 - prod(1-ep(USCaucs))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("exclusionPower", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fbiCaucs")
### * fbiCaucs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fbiCaucs
### Title: CODIS STR Loci allele frequency data
### Aliases: fbiCaucs
### Keywords: datasets

### ** Examples


data(fbiCaucs)
names(fbiCaucs)
fbiCaucs$loci
names(fbiCaucs$freqs)
fbiCaucs$freqs[[1]]
names(fbiCaucs$freqs[[1]])
fbiCaucs$freqs[[1]][1]




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fbiCaucs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fetchBMdata")
### * fetchBMdata

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fetchBMdata
### Title: Retrieve data from Budowle and Moretti (1999) from the web
### Aliases: fetchBMdata
### Keywords: datasets

### ** Examples


## not run
## Not run: 
##D db = fetchBMdata()
##D names(db)
##D f = db[["TRINIDADIAN"]]$freqs
##D dbExpect(f, k = "UN", collapse = TRUE)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fetchBMdata", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("locusIBS")
### * locusIBS

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: locusIBS
### Title: Identity by state at a locus
### Aliases: locusIBS

### ** Examples


data(fbiCaucs)
G = randomSample(1, fbiCaucs, rel = 'FS', N = 1000)
ibs = locusIBS(G)
barplot(tabulate(ibs+1, nbins = 3))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("locusIBS", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lrMix")
### * lrMix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lrMix
### Title: Calculate locuswise likelihood ratios for two person
###   victim/suspect mixtures
### Aliases: lrMix

### ** Examples


data(USCaucs)
p = randomProfilePairs(USCaucs, 10000)
log.lrs = log10(lrMix(p, USCaucs))
boxplot(log.lrs, las = 2)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lrMix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lrPC")
### * lrPC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lrPC
### Title: Likelihood Ratio for Parent-Child / Paternity Index
### Aliases: lrPC

### ** Examples


data(fbiCaucs)
P1 = randomProfile(fbiCaucs)
C1 = randomChild(P1, fbiCaucs)
lrPC(P1, C1, fbiCaucs)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lrPC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lrSib")
### * lrSib

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lrSib
### Title: Likelihood Ratio / Kinship Index for full-siblings
### Aliases: lrSib

### ** Examples


data(fbiCaucs)
P1 = randomProfile(fbiCaucs)
S1 = randomSib(P1, fbiCaucs)
P2 = randomProfile(fbiCaucs)
lrSib(P1, S1, fbiCaucs)
lrSib(P1, P2, fbiCaucs)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lrSib", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lrSibDebug")
### * lrSibDebug

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lrSibDebug
### Title: Likelihood Ratio / Kinship Index for full-siblings
### Aliases: lrSibDebug

### ** Examples


data(fbiCaucs)
P1 = randomProfile(fbiCaucs)
S1 = randomSib(P1, fbiCaucs)
P2 = randomProfile(fbiCaucs)
cat(paste(lrSibDebug(P1, S1, fbiCaucs)$Lines))
cat(paste(lrSibDebug(P1, P2, fbiCaucs)$Lines))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lrSibDebug", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("normalizeFreqs")
### * normalizeFreqs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: normalizeFreqs
### Title: Normalize frequencies to 1
### Aliases: normalizeFreqs

### ** Examples


data(fbiCaucs)

## induce an error
fbiCaucs$freqs[[1]] = rgamma(10,1,1)
checkFreqs(fbiCaucs)
fbiCaucs = normalizeFreqs(fbiCaucs)
checkFreqs(fbiCaucs)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("normalizeFreqs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.population")
### * print.population

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.population
### Title: Print summary details of a substructed population
### Aliases: print.population

### ** Examples


data(fbiCaucs)
p = breedFst(fbiCaucs)
print(p)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.population", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.profile")
### * print.profile

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.profile
### Title: Print a DNA profile
### Aliases: print.profile

### ** Examples


data(fbiCaucs)
P1 = randomProfile(fbiCaucs)
P2 = randomProfile(fbiCaucs)
P1
print(P1, horizontal = TRUE)
print(P2, horizontal = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.profile", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("randomChild")
### * randomChild

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: randomChild
### Title: Generate a random child from a given DNA profile and a given set
###   of allele frequencies
### Aliases: randomChild

### ** Examples


data(fbiCaucs)
P1 = randomProfile(fbiCaucs)
C1 = randomChild(P1,fbiCaucs)
P1
C1




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("randomChild", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("randomPCPairs")
### * randomPCPairs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: randomPCPairs
### Title: Generate one or more random parent/child pairs from a given set
###   of allele frequencies
### Aliases: randomPCPairs

### ** Examples


data(fbiCaucs)
P = randomPCPairs(fbiCaucs)
P$parent
P$child




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("randomPCPairs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("randomProfile")
### * randomProfile

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: randomProfile
### Title: Generate a random DNA profile from a given set of allele
###   frequencies
### Aliases: randomProfile

### ** Examples


data(fbiCaucs)
P1 = randomProfile(fbiCaucs)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("randomProfile", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("randomProfilePairs")
### * randomProfilePairs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: randomProfilePairs
### Title: Generate one or more random DNA profile pairs from a given set
###   of allele frequencies
### Aliases: randomProfilePairs

### ** Examples


data(fbiCaucs)
P = randomProfilePairs(fbiCaucs)
P$prof1
P$prof2




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("randomProfilePairs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("randomSample")
### * randomSample

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: randomSample
### Title: Generate a random sample of related (or unrelated) pairs of
###   people
### Aliases: randomSample

### ** Examples


data(fbiCaucs)
G = randomSample(1, fbiCaucs, "FS", 100)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("randomSample", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("randomSib")
### * randomSib

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: randomSib
### Title: Generate a random sibling from a given DNA profile and a given
###   set of allele frequencies
### Aliases: randomSib

### ** Examples


data(fbiCaucs)
P1 = randomProfile(fbiCaucs)
S1 = randomSib(P1,fbiCaucs)
P1
S1




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("randomSib", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("randomSibPairs")
### * randomSibPairs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: randomSibPairs
### Title: Generate one or more pairs of random siblings from a given set
###   of allele frequencies
### Aliases: randomSibPairs

### ** Examples


data(fbiCaucs)
P = randomSibPairs(fbiCaucs)
P$sib1
P$sib2




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("randomSibPairs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("readResults")
### * readResults

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: readResults
### Title: Read a simulation result set from file
### Aliases: readResults

### ** Examples


data(fbiCaucs)
## not run
## write the results of 100 unrelated profile pairs to
## results-sim-UN-100.csv.gz
## and read it back in
## Not run: 
##D sim(100, save = T)
##D unrel = readResults(100)
##D sim(100, rel = "FS", strVer = "01", save = T)
##D sibs = readResults(100, rel = "FS", strVer = "01")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("readResults", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sim")
### * sim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sim
### Title: Perform the relatives simulation
### Aliases: sim

### ** Examples


## not run
## this replicates Ge et al.'s experiment and takes about 45 minutes
## to run (I think)
## Not run: 
##D data(fbiCaucs)
##D N = 1000000
##D sim(N, fbiCaucs, save = T)
##D sim(N, fbiCaucs, 'FS', save = T)
##D sim(N, fbiCaucs, 'PC', save = T)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("toNexus")
### * toNexus

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: toNexus
### Title: Export a population with substructure to a Nexus file
### Aliases: toNexus

### ** Examples


data(USCaucs)
p = breedFst(USCaucs)
toNexus(p)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("toNexus", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("writeCSV")
### * writeCSV

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: writeCSV
### Title: Saves/writes population frequencies to disk
### Aliases: writeCSV

### ** Examples


data(USCaucs)
## Not run: 
##D   writeCSV("USCaucs.csv", USCaucs)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("writeCSV", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("writePop")
### * writePop

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: writePop
### Title: Saves/writes population profiles to disk
### Aliases: writePop

### ** Examples


data(USCaucs)
pop = breedFst(USCaucs)
## Not run: 
##D   writePop("USCaucs.csv", pop)
##D   
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("writePop", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
