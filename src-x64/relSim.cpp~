extern "C" {
  void locusLRmix(int *pnProfVic, int *pnProfSus, double *pdFreq, 
		  double *pdResult){ 
                  // double *pdF, int *nCase, double *pdResult){
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

  void profIbs(int *pnProf, int *nResult){
    *nResult = 0;

    int nA1 = pnProf[0];
    int nA2 = pnProf[1];
    int nB1 = pnProf[2];
    int nB2 = pnProf[3];

    if(nA1 == nB1 && nA2 == nB2){
      *nResult = 2;
    }else if((nA1 == nB1) || (nA2 == nB1) || (nA1 == nB2) || (nA2 == nB2)){
      *nResult = 1;
    } // else *nResult is zero
  }

  void LRmix(int *pnProfVic, int *pnProfSus, int *nLoci, double *pdFreqs,
	     int *pnOffsets, double *pdLocusLRs){
    int i;

    for(i = 0; i < *nLoci; i++){
      locusLRmix(&pnProfVic[2*i], &pnProfSus[2*i], &pdFreqs[pnOffsets[i]],
                 &pdLocusLRs[i]);
    }
  }

  void locusIbs(int *pnProfMat, int *pnResult, int *pN){

    // assumes pnProfMat is a vector of length 4*pN

    int i;
    for(i = 0; i < *pN; i++){
      int i1 = 4*i;
      int r;
      profIbs(&pnProfMat[i1], &r);
      pnResult[i] = r;
    }
  }

  void IBS(int *pnProf1, int *pnProf2, int *nLoci, int *nResult){
    int s = 0;
    int i;
    int *p = new int[4];

    for(i = 0; i < *nLoci; i++){
      int m;
      int i1 = 2*i;
      p[0] = pnProf1[i1];
      p[1] = pnProf1[i1+1];
      p[2] = pnProf2[i1];
      p[3] = pnProf2[i1+1];

      profIbs(p, &m);
      s +=  m;
    }

    delete [] p;

    *nResult = s;
  }

  void locusLRSib(int *pnProfSib1, int *pnProfSib2, double *pdFreq, double *pdResult){
    double dLR = 0;

    if(pnProfSib2[0]==pnProfSib2[1]){ // Sib2 is aa
      double dPa = pdFreq[pnProfSib2[0] - 1];

      if(pnProfSib1[0] == pnProfSib1[1]){ // Sib1 is aa or bb
	if(pnProfSib1[0] == pnProfSib2[0]){ // Sib1 is aa
	  dLR = (1+dPa)*(1+dPa)/(4*dPa*dPa);
	}else{ // Sib1 is bb
	  dLR = 0.25;
	}
      }else{ // Sib1 is ab or bc
	if(pnProfSib1[0]!=pnProfSib2[0] && pnProfSib1[1]!=pnProfSib2[0]){ // Sib1 is bc
	  dLR = 0.25;
	}else{ // Sib1 is ab
	  dLR = (1+dPa)/(4*dPa);
	}
      }
    }else{ // Sib2 is ab
      double dPa = pdFreq[pnProfSib2[0] - 1];
      double dPb = pdFreq[pnProfSib2[1] - 1];

      if(pnProfSib1[0]==pnProfSib1[1]){ // Sib1 is aa or bb or cc
	if(pnProfSib1[0]==pnProfSib2[0]){ // Sib1 is aa
	  dLR = (1+dPa)/(4*dPa);
	}else if(pnProfSib1[0]==pnProfSib2[1]){ // Sib1 is bb
	  dLR = (1+dPb)/(4*dPb);
	}else{ // Sib1 is cc
	  dLR = 0.25;
	}
      }else{ // Sib1 is ab, bc, ac, or cd
	if(pnProfSib1[0]!=pnProfSib2[0] && pnProfSib1[1]!=pnProfSib2[1] &&
	   pnProfSib1[1]!=pnProfSib2[0] && pnProfSib1[0]!=pnProfSib2[1]){ // Sib1 is cd
	  dLR = 0.25;
	}else if(pnProfSib1[0]==pnProfSib2[0] && pnProfSib1[1]==pnProfSib2[1]){ // Sib1 is ab
	  dLR = (1+dPa+dPb+2*dPa*dPb)/(8*dPa*dPb);
	}else{ // Sib1 is ac or bc
	  if((pnProfSib1[0]==pnProfSib2[0] && pnProfSib1[1]!=pnProfSib2[1]) ||
	     (pnProfSib1[1]==pnProfSib2[0] && pnProfSib1[0]!=pnProfSib2[1])){ // Sib1 is ac
	    dLR = (1+2*dPa)/(8*dPa);
	  }else{
	    dLR = (1+2*dPb)/(8*dPb);
	  }
	}
      }
    }

    *pdResult = dLR;
  }

  void lrSib(int *pnProfSib1, int *pnProfSib2, int *nLoci,
             double *pdFreq, int *pnAlleles, double *pdResult){

    double dLR = 1;
    int nOffset = 0;
    int nLoc;

    for(nLoc = 0; nLoc < *nLoci; nLoc++){
      double llr;
      int i1 = 2*nLoc;

      locusLRSib(&pnProfSib1[i1], &pnProfSib2[i1], &pdFreq[nOffset], &llr);

      dLR *= llr;
      nOffset += pnAlleles[nLoc];
    }

    *pdResult = dLR;
  }

  void locusLRPC(int* pnProfParent, int *pnProfChild,
                 double *pdFreq, double *pdResult){
    double dLR = 1;

    if(pnProfChild[0]==pnProfChild[1]){ // Child is aa
      double dPa = pdFreq[pnProfChild[0] - 1];

      if(pnProfParent[0]==pnProfParent[1]){ // Parent is aa or bb
	if(pnProfParent[0] == pnProfChild[0]){ // Parent is aa
	  dLR/=dPa;
	}else{ // Parent is bb
	  dLR = 0;
	}
      }else{ // Parent is ab or bc
	if(pnProfParent[0]!=pnProfChild[0] && pnProfParent[1]!=pnProfChild[0]){ // Parent is bc
	  dLR = 0;
	}else{ // Parent is ab
	  dLR/=2*dPa;
	}
      }
    }else{ // Child is ab
      double dPa = pdFreq[pnProfChild[0] - 1];
      double dPb = pdFreq[pnProfChild[1] - 1];

      if(pnProfParent[0]==pnProfParent[1]){ // Parent is aa or bb or cc
	if(pnProfParent[0]==pnProfChild[0]){ // Parent is aa
	  dLR/=2*dPa;
	}else if(pnProfParent[0]==pnProfChild[1]){ // Parent is bb
	  dLR/=2*dPb;
	}else{ // Parent is cc
	  dLR = 0;
	}
      }else{ // Parent is ab, bc, ac, or cd
	if(pnProfParent[0]==pnProfChild[0] && pnProfParent[1]==pnProfChild[1]){ // Parent is ab
	  dLR = (dPa+dPb)/(4*dPa*dPb);
	}else if(pnProfParent[0]==pnProfChild[0] || pnProfParent[1]==pnProfChild[0]){ // Parent is ac or ca
	  dLR/=4*dPa;
	}else if(pnProfParent[0]==pnProfChild[1] || pnProfParent[1]==pnProfChild[1]){ // Parent is bc or cb
	  dLR/=4*dPb;
	}else{ // Parent is cd
	  dLR = 0;
	}
      }
    }

    *pdResult = dLR;
  }

  void lrPC(int* pnProfParent, int *pnProfChild, int *nLoci,
            double *pdFreq, int *pnAlleles, double *pdResult){

    int nLoc = 0;
    double dLR = 1;
    int nOffset = 0;

    while(dLR > 0 && nLoc < *nLoci){
      double llr;
      int i1 = 2*nLoc;

      locusLRPC(&pnProfParent[i1], &pnProfChild[i1], &pdFreq[nOffset], &llr);

      dLR *= llr;
      nOffset += pnAlleles[nLoc];
      nLoc++;
    }

    *pdResult = dLR;

  }

  void locusProb(int *pnProf,double *pdFreq, double *pdResult){
    if(pnProf[0] == pnProf[1]){
      double dPa = pdFreq[pnProf[0] - 1];
      *pdResult = dPa*dPa;
    }else{
      double dPa = pdFreq[pnProf[0] - 1];
      double dPb = pdFreq[pnProf[1] - 1];
      *pdResult = 2*dPa*dPb;
    }
  }

  void prob(int* pnProf, int *nLoci,
            double *pdFreq, int *pnAlleles, double *pdResult){

    int nLoc;
    double dProd = 1;
    int nOffset = 0;

    for(nLoc = 0; nLoc < *nLoci; nLoc++){
      double lprob;
      int i1 = 2*nLoc;

      locusProb(&pnProf[i1], &pdFreq[nOffset], &lprob);

      dProd *= lprob;
      nOffset += pnAlleles[nLoc];
      nLoc++;
    }

    *pdResult = dProd;
  }
      
  int randAllele(double *pdFreqs, double dU){
    double dSum = pdFreqs[0];
    int nA = 0;

    while(dU > dSum){
      nA++;
      dSum += pdFreqs[nA];
    }

    return nA + 1;
  }

  void randomProfiles(int *pnProfile, int *nLoci, double *pdFreqs,
                      int *pnAlleles, double *pdU, int *nBlockSize){

    int nLoc;
    int i = 0;
    int nB;

    for(nB = 0; nB < *nBlockSize; nB++){
      int nOffset = 0;
      int nBlockOffset = nB*(*nLoci)*2;

      for(nLoc = 0; nLoc < *nLoci; nLoc++){
	int a1 = randAllele(&pdFreqs[nOffset], pdU[i++]);
	int a2 = randAllele(&pdFreqs[nOffset], pdU[i++]);

	if(a1 > a2){
	  pnProfile[nBlockOffset + 2*nLoc] = a2;
	  pnProfile[nBlockOffset + 2*nLoc + 1] = a1;
	}else{
	  pnProfile[nBlockOffset + 2*nLoc] = a1;
	  pnProfile[nBlockOffset + 2*nLoc + 1] = a2;
	}

	nOffset += pnAlleles[nLoc];
      }
    }
  }

  void randomSibs(int *pnProfSib1, int *pnProfSib2, int *nLoci,
                  double *pdFreqs, int *pnAlleles, double *pdU,
                  int *nBlockSize){

    int nLoc;
    int i = 0;
    int nB;

    for(nB = 0; nB < *nBlockSize; nB++){
      int nOffset = 0;
      int nBlockOffset = nB*(*nLoci)*2;

      for(nLoc = 0; nLoc < *nLoci; nLoc++){
	double dU = pdU[i++];
	int a1 = 0;
	int a2 = 0;

	if(dU < 0.25){
	  a1 = pnProfSib1[nBlockOffset + 2*nLoc];
	  a2 = pnProfSib1[nBlockOffset + 2*nLoc + 1];
	}else if(dU >= 0.25 && dU < 0.5){
	  a1 = pnProfSib1[nBlockOffset + 2*nLoc];
	  a2 = randAllele(&pdFreqs[nOffset], pdU[i++]);
	}else if(dU >= 0.50 && dU <= 0.75){
	  a1 = randAllele(&pdFreqs[nOffset], pdU[i++]);
  	  a2 = pnProfSib1[nBlockOffset + 2*nLoc + 1];
	}else{
	  a1 = randAllele(&pdFreqs[nOffset], pdU[i++]);
	  a2 = randAllele(&pdFreqs[nOffset], pdU[i++]);
	}

	if(a1 > a2){
	  pnProfSib2[nBlockOffset + 2*nLoc] = a2;
	  pnProfSib2[nBlockOffset + 2*nLoc + 1] = a1;
	}else{
	  pnProfSib2[nBlockOffset + 2*nLoc] = a1;
	  pnProfSib2[nBlockOffset + 2*nLoc + 1] = a2;
	}

	nOffset += pnAlleles[nLoc];
      }
    }
  }

  void randomChildren(int *pnProfParent, int *pnProfChild, int *nLoci,
                  double *pdFreqs, int *pnAlleles, double *pdU,
                  int *nBlockSize){

    int nLoc;
    int i = 0;
    int nB;

    for(nB = 0; nB < *nBlockSize; nB++){
      int nOffset = 0;
      int nBlockOffset = nB*(*nLoci)*2;

      for(nLoc = 0; nLoc < *nLoci; nLoc++){
	double dU = pdU[i++];
	int a1 = 0;
	int a2 = 0;

	if(dU < 0.5){
	  a1 = pnProfParent[nBlockOffset + 2*nLoc];
	  a2 = randAllele(&pdFreqs[nOffset], pdU[i++]);
	}else{
	  a1 = randAllele(&pdFreqs[nOffset], pdU[i++]);
	  a2 = pnProfParent[nBlockOffset + 2*nLoc + 1];
	}

	if(a1 > a2){
	  pnProfChild[nBlockOffset + 2*nLoc] = a2;
	  pnProfChild[nBlockOffset + 2*nLoc + 1] = a1;
	}else{
	  pnProfChild[nBlockOffset + 2*nLoc] = a1;
	  pnProfChild[nBlockOffset + 2*nLoc + 1] = a2;
	}

	nOffset += pnAlleles[nLoc];
      }
    }
  }

  void blockStatCounts(int *pnProf1, int *pnProf2, int *nLoci, int *nProf,
		       double *pdFreq, int *pnAlleles,
		       int *nCode, int *nFN,
                       int *pnIBSthresh, double *pdLRthresh,
		       int *pnResult, int *nNumResults){

    bool bFalseNeg = (bool)*nFN;
    int i, j;

    if(*nCode == 1){ // lrSib only
      for(i = 0; i < (*nProf); i++){
	int nOffset = 2*(*nLoci)*i;
	int *pProf1 = &pnProf1[nOffset];
	int *pProf2 = &pnProf2[nOffset];
	double dLRSib;

	lrSib(pProf1, pProf2, nLoci, pdFreq, pnAlleles, &dLRSib);

	if(bFalseNeg){
	  for(j = 0; j < *nNumResults; j++){
	    if(dLRSib < pdLRthresh[j])
	      pnResult[j] += 1;
	  }
	}else{
	  for(j = 0; j < *nNumResults; j++){
	    if(dLRSib >= pdLRthresh[j])
	      pnResult[j] += 1;
	  }
	}
      }
    }else if(*nCode == 2){ //lrPC only
      for(i = 0; i < (*nProf); i++){
	int nOffset = 2*(*nLoci)*i;
	int *pProf1 = &pnProf1[nOffset];
	int *pProf2 = &pnProf2[nOffset];
	double dLRPC;

	lrPC(pProf1, pProf2, nLoci, pdFreq, pnAlleles, &dLRPC);

	if(bFalseNeg){
	  for(j = 0; j < *nNumResults; j++){
	    if(dLRPC < pdLRthresh[j])
	      pnResult[j] += 1;
	  }
	}else{
	  for(j = 0; j < *nNumResults; j++){
	    if(dLRPC >= pdLRthresh[j])
	      pnResult[j] += 1;
	  }
	}
      }
    }else if(*nCode == 3){ // ibs only
      for(i = 0; i < (*nProf); i++){
	int nOffset = 2*(*nLoci)*i;
	int *pProf1 = &pnProf1[nOffset];
	int *pProf2 = &pnProf2[nOffset];
	int nIBS;

	IBS(pProf1, pProf2, nLoci, &nIBS);

	if(bFalseNeg){
	  for(j = 0; j < *nNumResults; j++){
	    if(nIBS < pnIBSthresh[j])
	      pnResult[j] += 1;
	  }
	}else{
	  for(j = 0; j < *nNumResults; j++){
	    if(nIBS >= pnIBSthresh[j])
	      pnResult[j] += 1;
	  }
	}
      }
    }else if(*nCode == 4){ // lrSib and ibs
      for(i = 0; i < (*nProf); i++){
	int nOffset = 2*(*nLoci)*i;
	int *pProf1 = &pnProf1[nOffset];
	int *pProf2 = &pnProf2[nOffset];
	double dLRSib;
        int nIBS;

	lrSib(pProf1, pProf2, nLoci, pdFreq, pnAlleles, &dLRSib);
	IBS(pProf1, pProf2, nLoci, &nIBS);

	if(bFalseNeg){
	  for(j = 0; j < *nNumResults; j++){
	    if(dLRSib < pdLRthresh[j] || nIBS < pnIBSthresh[j])
	      pnResult[j] += 1;
	  }
	}else{
	  for(j = 0; j < *nNumResults; j++){
	    if(dLRSib >= pdLRthresh[j] && nIBS >= pnIBSthresh[j])
	      pnResult[j] += 1;
	  }
	}
      }
    }else if(*nCode == 5){ // lrPC and ibs
      for(i = 0; i < (*nProf); i++){
	int nOffset = 2*(*nLoci)*i;
	int *pProf1 = &pnProf1[nOffset];
	int *pProf2 = &pnProf2[nOffset];
	double dLRPC;
        int nIBS;

	lrPC(pProf1, pProf2, nLoci, pdFreq, pnAlleles, &dLRPC);
	IBS(pProf1, pProf2, nLoci, &nIBS);

	if(bFalseNeg){
	  for(j = 0; j < *nNumResults; j++){
	    if(dLRPC < pdLRthresh[j] || nIBS < pnIBSthresh[j])
	      pnResult[j] += 1;
	  }
	}else{
	   for(j = 0; j < *nNumResults; j++){
	    if(dLRPC >= pdLRthresh[j] && nIBS >= pnIBSthresh[j])
	      pnResult[j] += 1;
	  }
	}
      }
    }else if(*nCode == 6){ // lrSib, lrPC and ibs

    }
  }

  void blockStats(int *pnProf1, int *pnProf2, int *nLoci, int *nProf,
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

	IBS(pProf1, pProf2, nLoci, &pnIBS[i]);
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
	IBS(pProf1, pProf2, nLoci, &pnIBS[i]);
      }
    }
  }

  void maximizeLRPC(int *pnProfParent, int *pnProfChild, int *nLoci,
                  double *pdFreqs, int *pnAlleles, double *pdU,
                  int *nBlockSize, double *dMax){

    int b;

    int *pnProf1 = new int[2*(*nLoci)*(*nBlockSize)];
    int *pnProf2 = new int[2*(*nLoci)*(*nBlockSize)];

    randomProfiles(pnProf1, nLoci, pdFreqs, pnAlleles, pdU, nBlockSize);

    int nOffset = 2*(*nLoci)*(*nBlockSize);
    randomChildren(pnProf1, pnProf2, nLoci, pdFreqs, pnAlleles, &pdU[nOffset],
                   nBlockSize);

    int iMax = 0;

    for(b = 0; b < *nBlockSize; b++){
      // calculate the LR

      double lr = 0;
      nOffset = 2*(*nLoci)*b;

      lrPC(&pnProf1[nOffset], &pnProf2[nOffset], nLoci, pdFreqs, pnAlleles, &lr);
      // update

      if(lr > *dMax){
	*dMax = lr;
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
}
