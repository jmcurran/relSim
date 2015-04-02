#include <Rcpp.h>
using namespace Rcpp;

int profIbsC(IntegerVector::const_iterator Prof){
  int nResult = 0;
  
  int nA1 = pnProf[0];
  int nA2 = pnProf[1];
  int nB1 = pnProf[2];
  int nB2 = pnProf[3];
  
  if(nA1 == nB1 && nA2 == nB2){
    *nResult = 2;
  }else if((nA1 == nB1) || (nA2 == nB1) || (nA1 == nB2) || (nA2 == nB2)){
    *nResult = 1;
  } // else *nResult is zero
  
  return nResult;
}

NumericVector LRmixC(IntegerVecotr ProfVic, IntegerVector ProfSus, List listFreqs){
  int i;
  int nLoci = listFreqs.size();
  NumericVector LocusLRs(nLoci);
  
  for(i = 0; i < *nLoci; i++){
    NumericVector Freqs = as<NumericVector>(listFreqs[i]);
    LocusLRs[i] = locusLRmixC(ProfVic[2 * i], ProfSus[2 * i], Freqs);
  }
  
  return LocusLRs;
}

IntegerVector locusIbsC(IntegerVector ProfMat, int N){
  
  // assumes pnProfMat is a vector of length 4 * N
  IntegerVector result(N);
  
  int i;
  for(i = 0; i < N; i++){
    int i1 = 4 * i;
    result[i] = profIbsC(ProfMat[i1]);
  }
  
  return result;
}

int IBSC(IntegerVector Prof1, IntegerVector pnProf2, int nLoci){
  int s = 0;
  int i;
  IntegerVector p(4);
  
  for(i = 0; i < nLoci; i++){
    int i1 = 2 * i;
    p[0] = Prof1[i1];
    p[1] = Prof1[i1 + 1];
    p[2] = Prof2[i1];
    p[3] = Prof2[i1 + 1];  
    s += profIbsC(p);
  }
  return s;
}
  
double locusLRmixC(IntegerVector::const_iterator ProfVic, IntegerVector::const_iterator ProfSus, 
                  NumericVector Freq){ 
  double dLR;
  double dPa, dPb, dPc, dPd;
  
  dPa = dPb = dPc =  dPd = 0;
  
  if(ProfVic[0] == ProfVic[1]){ // hom vic aa
    if(ProfSus[0] == ProfSus[1]){ // hom sus aa or bb
      if(ProfSus[0] == ProfVic[0]){ // sus aa
        // *nCase = 1;
        dPa = Freq[ProfVic[0] - 1];
        dLR = 1 / (dPa * dPa);
      }else{ // sus bb
        // *nCase = 2;
        dPa = Freq[ProfVic[0] - 1];
        dPb = Freq[ProfSus[0] - 1];
        dLR = 1 / (dPb * (2 * dPa + dPb));
      }
    }else{ // suspect ab, bc
      if(ProfSus[0] == ProfVic[0] || ProfSus[1] == ProfVic[0]){ // sus ab
        // *nCase = 4;
        dPa =  Freq[ProfVic[0] - 1];
        dPb =  Freq[ProfSus[ProfSus[0] == ProfVic[0] ? 1 : 0] - 1];
        dLR = 1 / (dPb * (2 * dPa + dPb));
      }else{ // suspect bc
        // *nCase = 5;
        dPb = Freq[ProfSus[0] - 1];
        dPc = Freq[ProfSus[1] - 1];
        dLR = 1/(2 * dPb * dPc);
      }
    }
  }else{ // vic ab
    if(ProfSus[0] == ProfSus[1]){ // hom sus aa or bb or cc
      if(ProfSus[0] == ProfVic[0] || ProfSus[0] == ProfVic[1]){ // sus aa or bb
        // *nCase = 5;
        dPa = Freq[ProfVic[0] - 1];
        dPb = Freq[ProfVic[1] - 1];
        double dp = dPa + dPb;
        dLR = 1 / (dp * dp);
      }else{ // suspect cc
        // *nCase = 6;
        dPa = Freq[ProfVic[0] - 1];
        dPb = Freq[ProfVic[1] - 1];
        dPc = Freq[ProfSus[0] - 1];
        dLR = 1 / (dPc * (2 * (dPa + dPb) + dPc));
      }
    }else{ // sus ab, ac or bc or cd
      if(ProfSus[0] == ProfVic[0] && ProfSus[1] == ProfVic[1]){ // sus ab
        // *nCase = 7;
        dPa = Freq[ProfSus[0] - 1];
        dPb = Freq[ProfSus[1] - 1];
        double dp = dPa + dPb;
        dLR = 1 / (dp * dp);
      }else if(ProfSus[0] == ProfVic[0] || ProfSus[1] == ProfVic[0]){ // sus ac 
        // *nCase = 8;
        dPa = Freq[ProfVic[0] - 1];
        dPb = Freq[ProfVic[1] - 1];
      
        if(ProfSus[0] == ProfVic[0]){
          dPc = Freq[ProfSus[1] - 1];
        }else{
          dPc = Freq[ProfSus[0] - 1];
        }
        dLR = 1 / (dPc * (2 * (dPa + dPb) + dPc));
      }else if(ProfSus[0] == ProfVic[1] || ProfSus[1] == ProfVic[1]){ // sus bc
        // *nCase = 9;
        dPa = Freq[ProfVic[0] - 1];
        if(ProfSus[0] == ProfVic[1]){
          dPb = Freq[ProfSus[0] - 1];
          dPc = Freq[ProfSus[1] - 1];
        }else{
          dPb = Freq[ProfSus[1] - 1];
          dPc = Freq[ProfSus[0] - 1];
        }
        dLR = 1 / (dPc * (2 * (dPa + dPb) + dPc));
      }else{ // sus cd
        // *nCase = 10;
        dPc = Freq[ProfSus[0] - 1];
        dPd = Freq[ProfSus[1] - 1];
        dLR = 1 / (2 * dPc * dPd);
      }
    }
  }
  //pdF[0] = dPa;
  //pdF[1] = dPb;
  //pdF[2] = dPc;
  //pdF[3] = dPd;
  return dLR;
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

List blockStats(int *pnProf1, int *pnProf2, int *nLoci, int *nProf,
double *pdFreq, int *pnAlleles, int *nCode,
int *pnIBS, double *pdLRSib, double *pdLRPC){
  
  int i;
  
  if(*nCode == 1){ // lrSib only
  for(i = 0; i < (*nProf); i++){
    int nOffset = 2*(*nLoci)*i;
    int *pProf1 = &pnProf1[nOffset];
    int *pProf2 = &pnProf2[nOffset];
    
    lrSib(pProf1, pProf2, nLoci, pdFreq, pnAlleles, &pdLRSib[i]);
  }
  }else if(*nCode == 2){ //lrPC only
  for(i = 0; i < (*nProf); i++){
    int nOffset = 2*(*nLoci)*i;
    int *pProf1 = &pnProf1[nOffset];
    int *pProf2 = &pnProf2[nOffset];
    
    lrPC(pProf1, pProf2, nLoci, pdFreq, pnAlleles, &pdLRPC[i]);
  }
  }else if(*nCode == 3){ // ibs only
  for(i = 0; i < (*nProf); i++){
    int nOffset = 2*(*nLoci)*i;
    int *pProf1 = &pnProf1[nOffset];
    int *pProf2 = &pnProf2[nOffset];
    
    IBSC(pProf1, pProf2, nLoci, &pnIBS[i]);
  }
  }else if(*nCode == 4){ // lrSib and lrPC
  for(i = 0; i < (*nProf); i++){
    int nOffset = 2*(*nLoci)*i;
    int *pProf1 = &pnProf1[nOffset];
    int *pProf2 = &pnProf2[nOffset];
    
    lrSib(pProf1, pProf2, nLoci, pdFreq, pnAlleles, &pdLRSib[i]);
    lrPC(pProf1, pProf2, nLoci, pdFreq, pnAlleles, &pdLRPC[i]);
    
  }
  }else if(*nCode == 5){ // lrSib, lrPC and ibs
  for(i = 0; i < (*nProf); i++){
    int nOffset = 2*(*nLoci)*i;
    int *pProf1 = &pnProf1[nOffset];
    int *pProf2 = &pnProf2[nOffset];
    
    lrSib(pProf1, pProf2, nLoci, pdFreq, pnAlleles, &pdLRSib[i]);
    lrPC(pProf1, pProf2, nLoci, pdFreq, pnAlleles, &pdLRPC[i]);
    IBSC(pProf1, pProf2, nLoci, &pnIBS[i]);
  }
  }
}
