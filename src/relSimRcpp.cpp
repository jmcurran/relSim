#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

double locusLRmix(IntegerVector::const_iterator ProfVic, IntegerVector::const_iterator ProfSus, 
                  NumericVector Freq){ 
  double dLR;
  double dPa, dPb, dPc, dPd;
  
  dPa = dPb = dPc =  dPd = 0;
  
  if(pnProfVic[0] == pnProfVic[1]){ // hom vic aa
    if(pnProfSus[0] == pnProfSus[1]){ // hom sus aa or bb
      if(pnProfSus[0] == pnProfVic[0]){ // sus aa
        // *nCase = 1;
        dPa = pdFreq[pnProfVic[0] - 1];
        dLR = 1 / (dPa * dPa);
      }else{ // sus bb
        // *nCase = 2;
        dPa = pdFreq[pnProfVic[0] - 1];
        dPb = pdFreq[pnProfSus[0] - 1];
        dLR = 1 / (dPb * (2 * dPa + dPb));
      }
    }else{ // suspect ab, bc
      if(pnProfSus[0] == pnProfVic[0] || pnProfSus[1] == pnProfVic[0]){ // sus ab
        // *nCase = 4;
        dPa =  pdFreq[pnProfVic[0] - 1];
        dPb =  pdFreq[pnProfSus[pnProfSus[0] == pnProfVic[0] ? 1 : 0] - 1];
        dLR = 1 / (dPb * (2 * dPa + dPb));
      }else{ // suspect bc
        // *nCase = 5;
        dPb = pdFreq[pnProfSus[0] - 1];
        dPc = pdFreq[pnProfSus[1] - 1];
        dLR = 1/(2 * dPb * dPc);
      }
    }
  }else{ // vic ab
    if(pnProfSus[0] == pnProfSus[1]){ // hom sus aa or bb or cc
      if(pnProfSus[0] == pnProfVic[0] || pnProfSus[0] == pnProfVic[1]){ // sus aa or bb
        // *nCase = 5;
        dPa = pdFreq[pnProfVic[0] - 1];
        dPb = pdFreq[pnProfVic[1] - 1];
        double dp = dPa + dPb;
        dLR = 1 / (dp * dp);
      }else{ // suspect cc
        // *nCase = 6;
        dPa = pdFreq[pnProfVic[0] - 1];
        dPb = pdFreq[pnProfVic[1] - 1];
        dPc = pdFreq[pnProfSus[0] - 1];
        dLR = 1 / (dPc * (2 * (dPa + dPb) + dPc));
      }
    }else{ // sus ab, ac or bc or cd
      if(pnProfSus[0] == pnProfVic[0] && pnProfSus[1] == pnProfVic[1]){ // sus ab
        // *nCase = 7;
        dPa = pdFreq[pnProfSus[0] - 1];
        dPb = pdFreq[pnProfSus[1] - 1];
        double dp = dPa + dPb;
        dLR = 1 / (dp * dp);
      }else if(pnProfSus[0] == pnProfVic[0] || pnProfSus[1] == pnProfVic[0]){ // sus ac 
        // *nCase = 8;
        dPa = pdFreq[pnProfVic[0] - 1];
        dPb = pdFreq[pnProfVic[1] - 1];
      
        if(pnProfSus[0] == pnProfVic[0]){
          dPc = pdFreq[pnProfSus[1] - 1];
        }else{
          dPc = pdFreq[pnProfSus[0] - 1];
        }
        dLR = 1 / (dPc * (2 * (dPa + dPb) + dPc));
      }else if(pnProfSus[0] == pnProfVic[1] || pnProfSus[1] == pnProfVic[1]){ // sus bc
        // *nCase = 9;
        dPa = pdFreq[pnProfVic[0] - 1];
        if(pnProfSus[0] == pnProfVic[1]){
          dPb = pdFreq[pnProfSus[0] - 1];
          dPc = pdFreq[pnProfSus[1] - 1];
        }else{
          dPb = pdFreq[pnProfSus[1] - 1];
          dPc = pdFreq[pnProfSus[0] - 1];
        }
        dLR = 1 / (dPc * (2 * (dPa + dPb) + dPc));
    }else{ // sus cd
      // *nCase = 10;
      dPc = pdFreq[pnProfSus[0] - 1];
      dPd = pdFreq[pnProfSus[1] - 1];
      dLR = 1 / (2 * dPc * dPd);
    }
  }
}
  //pdF[0] = dPa;
  //pdF[1] = dPb;
  //pdF[2] = dPc;
  //pdF[3] = dPd;
  *pdResult = dLR;
}
  

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

double locusLRPCC(IntegerVector::const_iterator ProfParent, IntegerVector::const_iterator ProfChild, NumericVector& Freq){
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

double lrPCC(IntegerVector::const_iterator ProfParent, IntegerVector::const_iterator ProfChild, 
            List listFreqs){
	int nLoci = listFreqs.size();
	int nLoc = 0;
	double dLR = 1;
	
	while(dLR > 0 && nLoc < nLoci){
		int i1 = 2 * nLoc;
		NumericVector Freqs = as<NumericVector>(listFreqs[nLoc]);	
    
		dLR *= locusLRPCC(ProfParent + i1, ProfChild + i1, Freqs);;
		nLoc++;
	}
	
	return(dLR);
}

// [[Rcpp::export]]
List maximizeLRPCC(List listFreqs, int nBlockSize){  
  int b;
  int nOffset;
  int nLoci = listFreqs.size();
  
  IntegerVector Prof1 = randomProfilesC(listFreqs, nBlockSize);
  IntegerVector Prof2 = randomChildrenC(Prof1, listFreqs, nBlockSize);
  int iMax = 0;
  double dMax = 0;
  
  for(b = 0; b < nBlockSize; b++){
    // calculate the LR
    
    nOffset = 2 * nLoci * b;
    double lr = lrPCC(Prof1.begin() + nOffset, Prof2.begin() + nOffset, listFreqs);
    // update
    
    if(lr > dMax){
      dMax = lr;
      iMax = b;
    }
  }
  
  nOffset = 2 *  nLoci * iMax;
  
  List listResult;
  
  listResult["parent"] = IntegerVector(Prof1.begin() + nOffset, Prof1.begin() + nOffset + 2 * nLoci);
  listResult["child"] = IntegerVector(Prof2.begin() + nOffset, Prof2.begin() + nOffset + 2 * nLoci);
  /*for(nLoc = 0; nLoc < nLoci; nLoc++){
    int i1 = 2*nLoc;
    pnProfParent[i1] = pnProf1[nOffset + i1];
    pnProfParent[i1+1] = pnProf1[nOffset + i1 + 1];
    pnProfChild[i1] = pnProf2[nOffset + i1];
    pnProfChild[i1+1] = pnProf2[nOffset + i1 + 1];
  }*/
  
  return listResult;
}
