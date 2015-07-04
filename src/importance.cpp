#include <Rcpp.h>
#include <map>
#include <string>
#include <sstream>
#include <vector>
using namespace Rcpp;
using namespace std;

class profileGenerator {
  vector<double> locusProbs;
  vector<double> cumProbs;
  int nAlleles;
  
  public:
  class Profile {
    protected:
    map<int, int> mapCounts;
    
    public:
    Profile(){};
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
    
    double prob(NumericVector& locusProbs, bool bLog = true){
      map<int, int>::iterator i = mapCounts.begin();
      double dSum = 0;
      
      while(i != mapCounts.end()){
        dSum += (i->second) * locusProbs[i->first];
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
  
  Profile randProf(int numContributors){
    Profile prof;
    NumericVector u = runif(2 * numContributors);
    
    for(int i = 0; i < 2 * numContributors; i++){
      int a = 1;
      while(u[i] >= cumProbs[a] && a < nAlleles){
       // Rprintf("%d %.7f %7f\n", a, u[i], cumProbs[a]);
        a++;
      }
      //Rprintf("%d %.7f %7f\n", a, u[i], cumProbs[a]);
      //Rprintf("Generated %d\n", a - 1);
      prof[a] += 1;
    }
    
    return(prof);
  }
  
  void printLocusProbs(){
    for(int i =  0; i < nAlleles; i++)
      Rprintf("%.7f %.7f\n", locusProbs[i], cumProbs[i]);
  }
};

NumericVector log(NumericVector& x){
  int nx = x.size();
  NumericVector r(nx);
  
  for(int i = 0; i < nx; i++){
    r[i] = log(x[i]);
  }
  
  return r;
}

// [[Rcpp::export(".importance")]]
double importance(NumericVector g, NumericVector g0, int numContributors, int nIterations) {
  profileGenerator gen(g0);
  
  //gen.printLocusProbs();
  
  double dSum = 0;
  NumericVector log_g = log(g);
  NumericVector log_g0 = log(g0);
  List listProfiles;
  
  for(int i = 0; i < nIterations; i++){
    profileGenerator::Profile p = gen.randProf(numContributors);
    //listProfiles.push_back(p.asList()); 
    if(p.numAlleles() <= 2)
      dSum +=exp(p.prob(log_g) - p.prob(log_g0));
  }

  return dSum / nIterations;
}

// [[Rcpp::export(".sampleWeights")]]
NumericVector sampleWeights(NumericVector g, NumericVector g0, int numContributors, int nIterations) {
  profileGenerator gen(g0);
  
  NumericVector log_g = log(g);
  NumericVector log_g0 = log(g0);
  NumericVector w(nIterations);
  List listProfiles;
  
  for(int i = 0; i < nIterations; i++){
    profileGenerator::Profile p = gen.randProf(numContributors);
    w[i] = exp(p.prob(log_g) - p.prob(log_g0));
  }

  return w;
}

// [[Rcpp::export(".tabulateN")]]
IntegerVector tabulateN(NumericVector g, int numContributors, int nIterations) {
  profileGenerator gen(g);
  
  IntegerVector n(2 * numContributors);
  List listProfiles;
  
  for(int i = 0; i < nIterations; i++){
    profileGenerator::Profile p = gen.randProf(numContributors);
    n[p.numAlleles()] += 1;
  }

  return n;
}