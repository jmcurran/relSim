#include <math.h>
#include <iostream>

extern "C" {
  int ix, iy, iz;

  unsigned long zunif(){
    ix = (171 * ix) % 30269;
    iy = (172 * iy) % 30307;
    iz = (170 * iz) % 30223;

    return ix + iy + iz;
  }

  double runif(){
    ix = (171 * ix) % 30269;
    iy = (172 * iy) % 30307;
    iz = (170 * iz) % 30223;

    double r = ix / 30269.0 + iy / 30307.0 + iz / 30323.0;
    return r - floor(r);
  }

  void breed(int *pnParents, int *ns, int *Ns, int *nLoci,
	     int *x, int *y, int *z){
    int s;
    int subpopStart;
    int i;
    int nNumAlleles = (*ns)*(*Ns)*2*(*nLoci);

    int *pnChildren = new int[nNumAlleles];

    ix = *x;
    iy = *y;
    iz = *z;

    for(s = 0; s < *ns; s++){
      subpopStart = s*(*Ns)*2*(*nLoci);

      //i2 = subpopStart + *Ns - 1;

      int p1, p2;
      int *pChild, *pParent1, *pParent2;

      int ctr = 0;
      unsigned int g, u = zunif();
      int loc;

      for(i = 0; i < (*Ns); i++){
	p1 = subpopStart + 2*(*nLoci)*floor(*Ns*runif());
	p2 = subpopStart + 2*(*nLoci)*floor(*Ns*runif());

	pChild = &pnChildren[subpopStart + (2*(*nLoci)*i)];
	pParent1 = &pnParents[p1];
	pParent2 = &pnParents[p2];
			     
	for(loc = 0; loc < *nLoci; loc++){
	  g = u & 1; // this pops a bit off u
	  u >>= 1; // and this shifts to the right so the bit is removed


	  // g will be 0 or 1
	  int a1  = pParent1[2*loc + g];
	  
	  g = u & 1; // this pops a bit off u
	  u >>= 1; // and this shifts to the right so the bit is removed
	  ctr += 2;

	  int a2  = pParent2[2*loc + g];
	  
	  if(a1 > a2){
	    int nSwap = a1;
	    a1 = a2;
	    a2 = nSwap;
	  }

	  pChild[2*loc] = a1;
	  pChild[2*loc + 1] = a2;

	  if(ctr == 16){
	    u = zunif();
	    ctr = 0;
	  }
	}
      }
    }

    for(i = 0; i < nNumAlleles; i++){
      pnParents[i] = pnChildren[i];
    }

    // preserve the random number seeds

    *x = ix;
    *y = iy;
    *z = iz;

    delete [] pnChildren;
  }
}
