extern "C" {
  void calculateAlleleFrequencies(int *pnPop, int *pnSubpopIdx,
                                  int N, int ns, int nLoci,
                                  int *pnNumLocusAlleles,
               		          double ***pppdAlleleFreqs,
                                  double ***pppdHom,
                                  int *pnSubPopSize){
    // Note - we can make this way more efficient
    // and if it's an issue we'll do it later.

    int r;
    int nLoc;


    for(r = 0; r < ns; r++){
      pnSubPopSize[r] = 0;

      // Clear the frequency array
      for(nLoc = 0; nLoc < nLoci; nLoc++){
	  int nAlleles = pnNumLocusAlleles[nLoc];
	  int nA;
	  for(nA = 0; nA < nAlleles; nA++){
	    pppdAlleleFreqs[r][nLoc][nA] = 0;
	    pppdHom[r][nLoc][nA] = 0;
	  }
      }
    }

    // iterate through the profiles subpop by subpop

    int i;

    for(i = 0; i < N; i++){
      r = pnSubpopIdx[i] - 1;
      pnSubPopSize[r]++;

      // tally up the allele counts
      int *pnProfile = &pnPop[i * 2 * nLoci];
      int nA1, nA2;
      int i1, i2;

      for(nLoc = 0; nLoc < nLoci; nLoc++){
	i1 = 2*nLoc;
	i2 = 2*nLoc + 1;

	nA1 = pnProfile[i1] - 1; // Alleles assumed to be 1..nA
	nA2 = pnProfile[i2] - 1;

	pppdAlleleFreqs[r][nLoc][nA1] += 1.0;
	pppdAlleleFreqs[r][nLoc][nA2] += 1.0;

	if(nA1 == nA2)
	  pppdHom[r][nLoc][nA1] += 1.0;

      }
    }

    for(r = 0; r < ns; r++){
      // divide each allele count by 2*nSubPopSize[r]
      for(nLoc = 0; nLoc < nLoci; nLoc++){
	int nA;
	for(nA = 0; nA < pnNumLocusAlleles[nLoc]; nA++){
	  pppdAlleleFreqs[r][nLoc][nA] /= 2.0*pnSubPopSize[r];
	  pppdHom[r][nLoc][nA] /= (double)pnSubPopSize[r];

	  //calculate and store average value for each allele count
	  //in pppdAlleleFreqs[ns][nLoc][nA]
	  //the last member of the pppdAlleleFreqs 3-d vector

	  double pi = (double)pnSubPopSize[r] / (double)N;
	  if(r==0)
	    pppdAlleleFreqs[ns][nLoc][nA] = pppdAlleleFreqs[r][nLoc][nA] * pi;
	  else
	    pppdAlleleFreqs[ns][nLoc][nA] += pppdAlleleFreqs[r][nLoc][nA] * pi;
	}
      }
    }
  }


  void calcTheta(int nLoci, int nSubPop, int *pnNumLocusAlleles,
                 int *pnSubPopSize,
                 double ***pppdAlleleFreqs, double ***pppdHom,
		 double *pdResult){

    int i;
    double dSum_ni = 0;
    double dSum_niSq = 0;

    for(i = 0; i < nSubPop; i++){
      dSum_ni += pnSubPopSize[i];
      dSum_niSq += pnSubPopSize[i] * pnSubPopSize[i];
    }

    double nc = (dSum_ni - (dSum_niSq / (double)dSum_ni)) / (nSubPop - 1);
    double nbar = (double)dSum_ni / (double)nSubPop;
    double dNumerator = 0;
    double dDenominator = 0;

    int nLoc;
    int nA;
    int r = nSubPop;

    for(nLoc = 0; nLoc < nLoci; nLoc++){
      double dLocusNumerator = 0;
      double dLocusDenominator = 0;

      int nAlleles = pnNumLocusAlleles[nLoc];

      for(nA = 0; nA < nAlleles; nA++){
	double dS1;
	double dS2,dS2a,dS2b,dS2c,dS2d,dS2e,dS2f,dS2g;
	double sASq=0;
	double pAbar = pppdAlleleFreqs[r][nLoc][nA];//0.5*pppdAlleleFreqs[r][nLoc][nA]/dSum_ni; // CHECK
	double HAbar = 0;
	double pAi;

	if(pAbar > 0){
	  for(i = 0;i < r;i++){
	    pAi = pppdAlleleFreqs[i][nLoc][nA];
	    sASq += pnSubPopSize[i] * (pAi - pAbar) * (pAi - pAbar);
	    HAbar += 2 * (pnSubPopSize[i])*(pAi - pppdHom[i][nLoc][nA]);
	  }

	  sASq /= (r - 1) * nbar;
	  HAbar /= dSum_ni;

	  dS1 = sASq - ((pAbar * (1 - pAbar) - sASq * (r - 1) / r -
                        0.25 * HAbar) / (nbar - 1));
	  dS2a = pAbar * (1 - pAbar);
	  dS2b = nbar / (r * (nbar - 1));
	  dS2c = r * (nbar - nc) / nbar;
	  dS2d = nbar - 1;
	  dS2e = (r - 1) * (nbar - nc);
	  dS2f = sASq / nbar;
	  dS2g = (nbar - nc) * HAbar / (4 * nc * nc);
	  dS2 =  dS2a - dS2b * (dS2c * dS2a - dS2f * (dS2d + dS2e) - dS2g);

	  dLocusNumerator += dS1;
	  dLocusDenominator += dS2;
	  pdResult[nLoc] = dLocusNumerator /  dLocusDenominator;

	  dNumerator += dS1;
	  dDenominator += dS2;

	}
      }
    }

    pdResult[nLoci] = dNumerator/dDenominator;
  }


  void calcFst(int *pnPop, int *pnSubPopIdx, int *N, int *ns,
               int *nLoci,
               int *pnNumLocusAlleles,
               double *pdResult){

    // Allocate memory for the allele frequencies
    // The dimension of this array is going to be (ns+1)*(nLoci)*(nAlleleMax)
    // way to create a 3-d array
    double ***pppdAlleleFreqs = new double**[*ns + 1];
    double ***pppdHom = new double**[*ns];
    int *pnSubPopSize = new int[*ns];

    int r, nLoc;

    for(r = 0; r <= *ns; r++){
      pppdAlleleFreqs[r] = new double*[*nLoci];
      if(r < *ns)
	pppdHom[r] = new double*[*nLoci];

      for(nLoc = 0; nLoc < *nLoci; nLoc++){
	int nAlleles = pnNumLocusAlleles[nLoc];
	pppdAlleleFreqs[r][nLoc] = new double[nAlleles];

	if(r < *ns)
	  pppdHom[r][nLoc] = new double[nAlleles];
      }
    }

    // calculate the allele frequencies

    calculateAlleleFrequencies(pnPop, pnSubPopIdx,
                               *N, *ns, *nLoci,
                               pnNumLocusAlleles,
               		       pppdAlleleFreqs,
                               pppdHom,
                               pnSubPopSize);

    // calcFst

    calcTheta(*nLoci, *ns, pnNumLocusAlleles,
                 pnSubPopSize,
                 pppdAlleleFreqs, pppdHom,
		 pdResult);


    // free the memory

    delete [] pnSubPopSize;

    for(r = 0; r <= *ns; r++){
       for(nLoc = 0; nLoc < *nLoci; nLoc++){
	 delete [] pppdAlleleFreqs[r][nLoc];

	if(r < *ns)
	  delete [] pppdHom[r][nLoc];
      }
       delete [] pppdAlleleFreqs[r];

       if(r < *ns)
	 delete [] pppdHom[r];
    }

    delete [] pppdAlleleFreqs;
    delete [] pppdHom;
  }

}
