extern "C" {
  
  
  
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
      
      IBSC(pProf1, pProf2, nLoci, &nIBS);
      
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
      IBSC(pProf1, pProf2, nLoci, &nIBS);
      
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
      IBSC(pProf1, pProf2, nLoci, &nIBS);
      
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
  

  

}
