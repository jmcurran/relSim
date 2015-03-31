#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
IntegerVector randomProfiles(List listFreqs, int nBlockSize){
  int nLoci = listFreqs.size();
  int nLoc;
  int nB;

  // pull in the sample function from R
  Environment base("package:base");
  Function sample = base["sample"];
  
  // reserve space for result
  IntegerVector Profiles(2 * nLoci * nBlockSize);
  
  for(nLoc = 0; nLoc < nLoci; nLoc++){
    NumericVector locusFreqs = as<NumericVector>(listFreqs[nLoc]);
    int nAlleles = locusFreqs.size();
    IntegerVector randAlleles = sample(seq_len(nAlleles), Named("size", 2 * nBlockSize),
                                       Named("replace", true),
                                       Named("prob", locusFreqs));
  
    for(nB = 0; nB < nBlockSize; nB++){
      int nBlockOffset = 2 * nLoci * nB;
    
      int a1 = randAlleles[2 * nB];
      int a2 = randAlleles[2 * nB + 1];
      
      if(a1 > a2){
        Profiles[nBlockOffset + 2 * nLoc] = a2;
        Profiles[nBlockOffset + 2 * nLoc + 1] = a1;
      }else{
        Profiles[nBlockOffset + 2 * nLoc] = a1;
        Profiles[nBlockOffset + 2 * nLoc + 1] = a2;
      }
    }
  }
  
  return Profiles;
}

IntegerVector randomSibs(IntegerVector ProfSib1, List listFreqs, int nBlockSize){
  int nLoci = listFreqs.size();
  int nLoc;
  int nB;
 
 // pull in the sample function from R
  Environment base("package:base");
  Function sample = base["sample"];
   
  // reserve space for the results
  IntegerVector ProfSib2(2 * nLoci * nBlockSize);
  
  
  // the mean + 10 sd rounded up should mean enough random alleles are available
  int sampleSize = (int)(100 * ceil(floor(nBlockSize + 10 * sqrt(nBlockSize * 0.5)) / 100));
    
 
  for(nLoc = 0; nLoc < nLoci; nLoc++){
    NumericVector U = runif(nBlockSize);
    
    NumericVector locusFreqs = as<NumericVector>(listFreqs[nLoc]);
    int nAlleles = locusFreqs.size();
    IntegerVector randAlleles = sample(seq_len(nAlleles), Named("size", sampleSize),
                                       Named("replace", true),
                                       Named("prob", locusFreqs));
    
    int j = 0;
    int a1, a2;
    
    for(nB = 0; nB < nBlockSize; nB++){
      int nBlockOffset = 2 * nLoci * nB;
    
      if(U[nB] < 0.25){
        a1 = ProfSib1[nBlockOffset + 2 * nLoc];
        a2 = ProfSib1[nBlockOffset + 2 * nLoc + 1];
      }else if(U[nB] >= 0.25 && U[nB] < 0.5){
        a1 = ProfSib1[nBlockOffset + 2*nLoc];
        a2 = randAlleles[j++];
      }else if(U[nB] >= 0.50 && U[nB] <= 0.75){
        a1 = randAlleles[j++];
        a2 = ProfSib1[nBlockOffset + 2 * nLoc + 1];
      }else{
        a1 = randAlleles[j++];
        a2 = randAlleles[j++];
      }
    
      if(a1 > a2){
        ProfSib2[nBlockOffset + 2*nLoc] = a2;
        ProfSib2[nBlockOffset + 2*nLoc + 1] = a1;
      }else{
        ProfSib2[nBlockOffset + 2*nLoc] = a1;
        ProfSib2[nBlockOffset + 2*nLoc + 1] = a2;
      }  
    }
  }
  
  return ProfSib2;
}

IntegerVector randomChildren(IntegerVector ProfParent, List listFreqs, int nBlockSize){
  int nLoci = listFreqs.size();
  int nLoc;
  int nB;
 
 // pull in the sample function from R
  Environment base("package:base");
  Function sample = base["sample"];
   
  // reserve space for the results
  IntegerVector ProfChild(2 * nLoci * nBlockSize);
  
  int sampleSize = (int)ceil(nBlockSize * 0.5);
  
   for(nLoc = 0; nLoc < nLoci; nLoc++){
    NumericVector U = runif(nBlockSize);
    
    NumericVector locusFreqs = as<NumericVector>(listFreqs[nLoc]);
    int nAlleles = locusFreqs.size();
    IntegerVector randAlleles = sample(seq_len(nAlleles), Named("size", sampleSize),
                                       Named("replace", true),
                                       Named("prob", locusFreqs));
    
    int j = 0;
    int a1, a2;
    
    for(nB = 0; nB < nBlockSize; nB++){
       int nBlockOffset = 2 * nLoci * nB;
     
      
      if(U[nB] < 0.5){
        a1 = ProfParent[nBlockOffset + 2 * nLoc];
        a2 = randAlleles[j++];
      }else{
        a1 = randAlleles[j++];
        a2 = ProfParent[nBlockOffset + 2 * nLoc + 1];
      }
      
      if(a1 > a2){
        ProfChild[nBlockOffset + 2*nLoc] = a2;
        ProfChild[nBlockOffset + 2*nLoc + 1] = a1;
      }else{
        ProfChild[nBlockOffset + 2*nLoc] = a1;
        ProfChild[nBlockOffset + 2*nLoc + 1] = a2;
      }
    }
  }
  
  return ProfChild;
}