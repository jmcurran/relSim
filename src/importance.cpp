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
    int numTotalAlleles;
    long numTotalAllelesFactorial;
    map<int, int> mapCounts;
    
    public:
    Profile(const profileGenerator& pg, int nTotalAlleles)
      : nAlleles(pg.nAlleles){
      numTotalAlleles = nTotalAlleles;
      numTotalAllelesFactorial = factorial(numTotalAlleles);
    };
      
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
        result[i->first - 1] = i->second;
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
      double dFact = std::log(numTotalAllelesFactorial);
      
      while(i != mapCounts.end()){
        dSum += (i->second) * std::log(locusProbs[i->first - 1]);
        dFact -= std::log(i->second);
        i++;
      }
      
      return dSum + dFact;
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
    Profile prof(*this, 2 * numContributors);
    
    // choose numAllelesShowing alleles without replacement
    IntegerVector alleles = seq_len(locusProbs.size());
    
    if(numAllelesShowing != locusProbs.size()){
      alleles = sample(alleles, numAllelesShowing, false, locusProbs);
    }
    
    for(IntegerVector::iterator a = alleles.begin(); a != alleles.end(); a++){
      prof[*a] = 1;
     // Rprintf("%d\n",*a);
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
List IS(NumericVector freqs,int N, int numContributors, int numAllelesShowing){
  profileGenerator g(freqs);
  
  NumericMatrix Alleles(N, freqs.size());
  NumericVector probs(N);
  
  for(int i = 0; i < N; i++){
    profileGenerator::Profile p = g.randProf(numContributors, numAllelesShowing);
    Alleles(i, _) = p.asNumericVector();
    probs[i] = p.prob(freqs);
  }
  
  List results;
  results["Alleles"] = Alleles;
  results["probs"] = probs;
  
  return results;
}

// [[Rcpp::export]]
NumericVector ISprob(const NumericVector& freqs, const NumericMatrix& AlleleCombs, const NumericMatrix& Perms){
  int numCombs = AlleleCombs.nrow();
  int numAlleles = AlleleCombs.ncol();
  int numPerms = Perms.nrow();
  NumericVector results(numCombs);
  
  for(int i = 0; i < numCombs; i++){
   
    for(int j = 0; j < numPerms; j++){
   
      double p = freqs[AlleleCombs(i, Perms(j,0) - 1) - 1];
      double s = p;
      
      for(int k = 1; k < numAlleles; k++){
        double pk = freqs[AlleleCombs(i, Perms(j, k) - 1) - 1];
        p *=  pk / (1 - s);
        s += pk;
      }
   
      results[i] += p;
    }
  }
  
  return results;
}
