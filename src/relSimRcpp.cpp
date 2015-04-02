#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
IntegerVector randomProfilesC(List listFreqs, int nBlockSize){
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

// [[Rcpp::export]]
IntegerVector randomSibsC(IntegerVector ProfSib1, List listFreqs, int nBlockSize){
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

// [[Rcpp::export]]
IntegerVector randomChildrenC(IntegerVector ProfParent, List listFreqs, int nBlockSize){
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
        ProfChild[nBlockOffset + 2 * nLoc] = a2;
        ProfChild[nBlockOffset + 2 * nLoc + 1] = a1;
      }else{
        ProfChild[nBlockOffset + 2 * nLoc] = a1;
        ProfChild[nBlockOffset + 2 * nLoc + 1] = a2;
      }
    }
  }
  
  return ProfChild;
}

double locusLRPCC(IntegerVector& ProfParent, IntegerVector& ProfChild, NumericVector& Freq, int nLoc){
	double dLR = 1;
	
	if(ProfChild[0] == ProfChild[1]){ // Child is aa
		double dPa = Freq[ProfChild[0] - 1];
  
		if(ProfParent[0] == ProfParent[1]){ // Parent is aa or bb
			if(ProfParent[0] == ProfChild[0]){ // Parent is aa
				dLR /= dPa;
			}else{ // Parent is bb
				dLR = 0;
			}
		}else{ // Parent is ab or bc
			if(ProfParent[0] != ProfChild[0] && ProfParent[1] != ProfChild[0]){ // Parent is bc
				dLR = 0;
			}else{ // Parent is ab
				dLR /= 2 * dPa;
			}
		}
	}else{ // Child is ab
		double dPa = Freq[ProfChild[0] - 1];
		double dPb = Freq[ProfChild[1] - 1];
  
		if(ProfParent[0]==ProfParent[1]){ // Parent is aa or bb or cc
			if(ProfParent[0]==ProfChild[0]){ // Parent is aa
				dLR /= 2 * dPa;
			}else if(ProfParent[0]==ProfChild[1]){ // Parent is bb
				dLR /= 2 * dPb;
			}else{ // Parent is cc
				dLR = 0;
			}
		}else{ // Parent is ab, bc, ac, or cd
			if(ProfParent[0] == ProfChild[0] && ProfParent[1] == ProfChild[1]){ // Parent is ab
				dLR = (dPa + dPb) / (4 * dPa * dPb);
			}else if(ProfParent[0] == ProfChild[0] || ProfParent[1] == ProfChild[0]){ // Parent is ac or ca
				dLR /= 4 * dPa;
			}else if(ProfParent[0] == ProfChild[1] || ProfParent[1] == ProfChild[1]){ // Parent is bc or cb
				dLR /= 4 * dPb;
			}else{ // Parent is cd
				dLR = 0;
			}
		}
	}
	
  return(dLR);
}

double lrPC(IntegerVector ProfParent, IntegerVector ProfChild, List listFreqs, int nPos){
	int nLoci = listFreqs.size();
	int nLoc = 0;
	double dLR = 1;
	int nOffset = 0;
	
	while(dLR > 0 && nLoc < nLoci){
		double llr;
		

		NumericVector Freqs = as<NumericVector>(listFreqs[nLoc]);
		
		locusLRPC(ProfParent, ProfChild, Freqs, nLoc, nPos);
		
		dLR *= llr;
		nLoc++;
	}
	
	return(dLR);
}

// [[Rcpp::export]]
double maximizeLRPC(IntegerVector ProfParent, IntegerVector ProfChild, 
                    List listFreqs, int nBlockSize){
  
  double dMax;
  int b;
  int nOffset;
  
  IntegerVector Prof1 = randomProfiles(listFreqs, nBlockSize);
  IntegerVector Prof2 = randomChildren(Prof1, listFreqs, nBlockSize);
  int iMax = 0;
  
  for(b = 0; b < nBlockSize; b++){
    // calculate the LR
    
    double lr = 0;
    nOffset = 2 * nLoci * b;
    
    lrPC(&pnProf1[nOffset], &pnProf2[nOffset], nLoci, pdFreqs, pnAlleles, &lr);
    // update
    
    if(lr > dMax){
      dMax = lr;
      iMax = b;
    }
  }
  
  int nLoc;
  nOffset = 2*(*nLoci)*iMax;
  
  for(nLoc = 0; nLoc < *nLoci; nLoc++){
    int i1 = 2*nLoc;
    pnProfParent[i1] = pnProf1[nOffset + i1];
    pnProfParent[i1+1] = pnProf1[nOffset + i1 + 1];
    pnProfChild[i1] = pnProf2[nOffset + i1];
    pnProfChild[i1+1] = pnProf2[nOffset + i1 + 1];
  }
  
  
  delete [] pnProf1;
  delete [] pnProf2;
}


void f1(vector<double>& x){
	cout << x.size() endl;
}

void f2(void){
	double z[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	vector<double> y(z, z + 10);
	
	for(int i = 0; i < 10; i++)
		f1(y+i);
	
}