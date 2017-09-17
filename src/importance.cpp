#include <Rcpp.h>
#include <math.h>
#include <map>
#include <string>
#include <sstream>
#include <vector>
using namespace Rcpp;
using namespace std;

class profileGenerator {
  NumericVector locusProbs;
  vector<double> cumProbs;
  int nAlleles;
  
  public:
  class Profile {
    protected:
    int nAlleles;
    map<int, int> mapCounts;
    
    public:
    Profile(const profileGenerator& pg)
      : nAlleles(pg.nAlleles){};
    int numAlleles(){
      return mapCounts.size();
    }
    
    // get
    const int operator[](int a) const{
      return mapCounts.at(a);
    };
    
    // set
    int& operator[](int a){
      return mapCounts[a];
    };
 
    List asList(){
      List result;
      map<int, int>::iterator i = mapCounts.begin();
      ostringstream oss;
     
      while(i != mapCounts.end()){
        oss << i->first;
        result[oss.str().c_str()] = i->second;
        oss.str("");
        i++;
      }
      
      return result;
    }
    
    NumericVector asNumericVector(){
      NumericVector result(nAlleles);
      map<int, int>::iterator i = mapCounts.begin();
      
      while(i != mapCounts.end()){
        result[i->first] = i->second;
        i++;
      }
      
      return result;
    }
    
    long factorial(long n){
      switch(n){
      case 0:
      case 1:
        return 1;
      case 2:
        return 2;
      case 3:
        return 6;
      case 4:
        return 24;
      case 5:
        return 120;
      case 6:
        return 720;
      case 7:
        return 5040;
      case 8:
        return 40320;
      case 9:
        return 362880;
      case 10:
        return 3628800;
      case 11:
        return 39916800;
      case 12:
        return 479001600;
      default:
        return n * factorial(n - 1);
      }
    }
    
    double prob(NumericVector& locusProbs, bool bLog = true){
      map<int, int>::iterator i = mapCounts.begin();
      double dSum = 0;
      
      while(i != mapCounts.end()){
        dSum += (i->second) * std::log(locusProbs[i->first]);
        i++;
      }
      
      return dSum;
    }
    
  };
  
  
  profileGenerator(NumericVector& locus){
    nAlleles = locus.size();
    
    for(int a = 0; a < nAlleles; a++){
      locusProbs.push_back(locus[a]);
      if(a == 0){
        cumProbs.push_back(locus[0]);
      }else{
        cumProbs.push_back(cumProbs[a - 1] + locus[a]);
      }
    }
  };
  
  Profile randProf(int numContributors, int numAllelesShowing){
    Profile prof(*this);
    
    // choose numAllelesShowing alleles without replacement
    IntegerVector alleles = sample((int)locusProbs.size(), numAllelesShowing, false, locusProbs, false);
    
    for(IntegerVector::iterator a = alleles.begin(); a != alleles.end(); a++){
      prof[*a] = 1;
    }
    
    alleles = sample(alleles, 2 * (numContributors) - numAllelesShowing, true);
    
    for(IntegerVector::iterator a = alleles.begin(); a != alleles.end(); a++){
      prof[*a] += 1;
    }
    
    return(prof);
  }
  
  void printLocusProbs(){
    for(int i =  0; i < nAlleles; i++)
      Rprintf("%.7f %.7f\n", locusProbs[i], cumProbs[i]);
  }
};

// [[Rcpp::export]]  
NumericMatrix IS(NumericVector freqs,int N, int numContributors, int numAllelesShowing){
  profileGenerator g(freqs);
  
  NumericMatrix result(N, freqs.size());
  
  for(int i = 0; i < N; i++){
    result(i, _) = g.randProf(numContributors, numAllelesShowing).asNumericVector();
  }
  
  return result;
}

