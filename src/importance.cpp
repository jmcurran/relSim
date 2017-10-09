#include <Rcpp.h>
#include <math.h>
#include <map>
#include <string>
#include <sstream>
#include <vector>
using namespace Rcpp;
using namespace std;

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

class Profile {
public:
  class Locus {
  protected:
    int nAlleles;
    int numTotalAlleles;
    long numTotalAllelesFactorial;
    map<int, int> mapCounts;

  public:
    Locus(int numAlleles, int nTotalAlleles){
      nAlleles = numAlleles;
      numTotalAlleles = nTotalAlleles;
      numTotalAllelesFactorial = factorial(numTotalAlleles);
    }
    
    Locus(const Locus& loc){
      nAlleles = loc.nAlleles;
      numTotalAlleles = loc.numTotalAlleles;
      numTotalAllelesFactorial = loc.numTotalAllelesFactorial;
      mapCounts = loc.mapCounts;
    }
    
    Locus& operator=(const Locus& loc){
      nAlleles = loc.nAlleles;
      numTotalAlleles = loc.numTotalAlleles;
      numTotalAllelesFactorial = loc.numTotalAllelesFactorial;
      mapCounts = loc.mapCounts;
      
      return *this;
    }

    int numAlleles(){
      return mapCounts.size();
    }

    // get
    const int operator[](int a) const{
      return mapCounts.at(a);
    }

    // set
    int& operator[](int a){
      return mapCounts[a];
    }

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

    double prob(NumericVector& locusProbs, bool bLog = true){
      map<int, int>::iterator i = mapCounts.begin();
      double dSum = 0;

      if(mapCounts.size() == 1)
        return (i->second) * std::log(locusProbs[i->first - 1]);

      double dFact = std::log(numTotalAllelesFactorial);

      while(i != mapCounts.end()){
        dSum += (i->second) * std::log(locusProbs[i->first - 1]);
        dFact -= std::log(factorial(i->second));
        i++;
      }

      return dSum + dFact;
    }

  //   Locus randLoc(int numContributors, int numAllelesShowing){
  //     Locus loc(*this, 2 * numContributors);
  // 
  //     // choose numAllelesShowing alleles without replacement
  //     IntegerVector alleles = seq_len(probs.size());
  // 
  //     if(numAllelesShowing != probs.size()){
  //       alleles = sample(alleles, numAllelesShowing, false, locusProbs);
  //     }
  // 
  //     double sum = 0;
  //     NumericVector newProbs;
  // 
  //     for(IntegerVector::iterator a = alleles.begin(); a != alleles.end(); a++){
  //       int allele = *a;
  //       prof[allele] = 1;
  //       double f = locusProbs[allele - 1];
  //       newProbs.push_back(f);
  //       sum += f;
  //       // Rprintf("%d\n",*a);
  //     }
  // 
  //     // normalize
  //     newProbs = newProbs / sum;
  // 
  //     alleles = sample(alleles, 2 * (numContributors) - numAllelesShowing, true, newProbs);
  // 
  //     for(IntegerVector::iterator a = alleles.begin(); a != alleles.end(); a++){
  //       prof[*a] += 1;
  //     }
  // 
  //     return(prof);
  //   }
  // 
  // };
  };

  vector<Locus> profile;
  
  Profile(int numLoci, int numContributors, const IntegerMatrix& numAllelesShowing){
    for(int loc = 0; loc < numLoci; loc++){
      profile.push_back(Locus(numAllelesShowing[loc], 2 * numContributors));
    }
  }
  
  Profile(const Profile& prof){
    profile = prof.profile;
  }
  
  Profile& operator=(const Profile& prof){
    profile = prof.profile;
    
    return *this;
  }
};
  

class profileGenerator {
public:
  class FreqInfo {
    vector<double> probs;
    vector<double> cumProbs;
    int numAlleles;
    
  public:
    FreqInfo(const vector<double>& lprobs){
      numAlleles = lprobs.size();
      
      vector<double>::const_iterator i = lprobs.begin();
      int a = 0;
      
      while(i != lprobs.end()){
        double f = *i;
        probs.push_back(f);
        
        if(a == 0){
          cumProbs.push_back(f);
        }else{
          cumProbs.push_back(cumProbs[a - 1] + f);
        }
        a++; 
        i++;
      }
    }
    
    FreqInfo& operator=(const vector<double>& lprobs){
      numAlleles = lprobs.size();
      
      vector<double>::const_iterator i = lprobs.begin();
      int a = 0;
      
      while(i != lprobs.end()){
        double f = *i;
        probs.push_back(f);
        
        if(a == 0){
          cumProbs.push_back(f);
        }else{
          cumProbs.push_back(cumProbs[a - 1] + f);
        }
        a++; 
        i++;
      }
      
      return *this;
    }
    
    void print(){
      for(int i =  0; i < numAlleles; i++){
        Rprintf("%d: %.7f %.7f\n", i + 1, probs[i], cumProbs[i]);
      }
    }
  };
  
  vector<FreqInfo> freqs;
  int numLoci;

  
  profileGenerator(const List& freqs){
    numLoci = (int)freqs.size();
    List::const_iterator i = freqs.begin();
    
    while(i != freqs.end()){
      vector<double> f = as< vector<double> >((NumericVector)*i);
      this->freqs.push_back(f);
      i++;
    }
  };
  
  void print(void){
    vector<FreqInfo>::iterator i = freqs.begin();
    int loc = 1;
    
    while(i != freqs.end()){
      Rprintf("Locus %2d\n", loc);
      Rprintf("-----------\n");
      i->print();
      Rprintf("-----------\n\n");
      i++;
      loc++;
    }
  }
  
  // [[Rcpp::plugins("cpp11")]]
  
  Profile randProf(int numContributors, const IntegerVector& numAllelesShowing){
    Profile prof(numLoci, numContributors, numAllelesShowing);
   
  }
  
  // void printLocusProbs(){
  //   for(int i =  0; i < nAlleles; i++)
  //     Rprintf("%.7f %.7f\n", locusProbs[i], cumProbs[i]);
  // }
};


// [[Rcpp::export(".IS")]]  
void IS(List freqs,int N, int numContributors, int maxAllelesShowing){
  profileGenerator g(freqs);
  int numLoci = g.numLoci;

  IntegerVector numAllelesShowing = sample(maxAllelesShowing, N * g.numLoci, TRUE);
  
  List profiles(N);
  IntegerVector::iterator nA = numAllelesShowing.begin();
  
  for(int i = 0; i < N; i++){
    Profile p = g.randProf(numContributors, *(nA + i * numLoci));
  }
  
  
  // NumericMatrix Alleles(N, freqs.size());
  // NumericVector probs(N);
  // IntegerVector numAllelesShowing(N);
  // 
  // numAllelesShowing = sample(maxAllelesShowing, N, TRUE);
  // 
  // for(int i = 0; i < N; i++){
  //   profileGenerator::Profile p = g.randProf(numContributors, numAllelesShowing[i]);
  //   Alleles(i, _) = p.asNumericVector();
  //   probs[i] = p.prob(freqs);
  // }
  // 
  // List results;
  // results["Alleles"] = Alleles;
  // results["probs"] = probs;
  // 
  // return results;
}

// // [[Rcpp::export]]
// NumericVector ISprob(const List& listCombs, const List& Perms){
//   int numCombs = listCombs.size();
//   NumericVector results(numCombs);
//   
//   for(int i = 0; i < numCombs; i++){
//     List lComb = as<List>(listCombs[i]);
//     NumericVector freqs = as<NumericVector>(lComb["f"]);
//     IntegerVector alleles = as<IntegerVector>(lComb["a"]);
//     IntegerVector counts = as<IntegerVector>(lComb["c"]);
//     int numAlleles = as<int>(lComb["n"]);
//     
//     NumericMatrix perms = as<NumericMatrix>(Perms[numAlleles - 1]);
//     int numPerms = perms.nrow();
//    
//     for(int j = 0; j < numPerms; j++){
//       double p , s;
//       p = s = freqs[perms(j, 0) - 1];
//       
//       for(int k = 1; k < numAlleles; k++){
//         double pk = freqs[perms(j, k) - 1];
//         p *=  pk / (1 - s);
//         s += pk;
//       }
//      // Rprintf("p = %.7f\n", p);
//       results[i] += p;
//     }
//     
//    // Rprintf("%.7f\n", results[i]);
//     freqs = freqs / sum(freqs);
// 
//     int sumCounts = 0;
//     double p2 = 0;
// 
//     for(int j = 0; j < numAlleles; j++){
//       // Rprintf("%d %d %.7f %.7f\n", alleles[j], counts[j], freqs[j], p2);
//       p2 += (counts[j] - 1) * std::log(freqs[j]) - std::log(factorial(counts[j] - 1));
//       sumCounts += counts[j] - 1;
//     }
// 
//     p2 += std::log(factorial(sumCounts));
//     // Rprintf("%.7f\n", p2);
//     results[i] = std::log(results[i]) + p2;
//   }
//   
//   return results;
// }
