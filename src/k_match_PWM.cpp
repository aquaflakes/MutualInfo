//assign kmers to best scoring PWM

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
// #include "common_op.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdio>


// template< class C >
// auto cbegin( C& c ) -> decltype(c.cbegin());
//
// template< class C >
// auto cbegin( const C& c ) -> decltype(c.cbegin());



using namespace Rcpp;
using namespace std;


inline char complement (const char c)
{
  switch (c) {
  case 'N':case 'n':
    return 'N';
  case 'A':case 'a':
    return 'T';
  case 'G':case 'g':
    return 'C';
  case 'C':case 'c':
    return 'G';
  case 'T':case 't':
    return 'A';
  default:
    stop("Invalid nucleotide.");
  }
}


// [[Rcpp::export]]
std::string revComp(std::string seq)
{
  std::transform(seq.cbegin(), seq.cend(), seq.begin(), complement);
  reverse(seq.begin(),seq.end());
  return(seq);
}

template <int RTYPE>
IntegerVector fast_factor_template( const Vector<RTYPE>& x ) {
  Vector<RTYPE> levs = sort_unique(x);
  IntegerVector out = match(x, levs);
  out.attr("levels") = as<CharacterVector>(levs);
  out.attr("class") = "factor";
  return out;
}


// [[Rcpp::export]]
SEXP fast_factor( SEXP x ) {
  switch( TYPEOF(x) ) {
  case INTSXP: return fast_factor_template<INTSXP>(x);
  case REALSXP: return fast_factor_template<REALSXP>(x);
  case STRSXP: return fast_factor_template<STRSXP>(x);
  }
  return R_NilValue;
}


inline int findRow(char c)
{
  switch (c) {
  case 'A':case 'a':
    return 0;
  case 'C':case 'c':
    return 1;
  case 'G':case 'g':
    return 2;
  case 'T':case 't':
    return 3;
  default:
    stop("Invalid nucleotide");
  }
}

inline void scoring_k_with_PWM(std::string &kmer, NumericMatrix& currMotif, double &MaxScore, std::string &MotifName, double threshold,  std::string &currMotifName)
  {
    int posNo= currMotif.ncol()-kmer.length()+1;
    for (int pos=0; pos< posNo; pos++)
      {
        double scoreCurrPos=0;
        for (int kpos=0; kpos< kmer.length(); kpos++)
          scoreCurrPos+= currMotif(findRow(kmer[kpos]),pos+kpos);
          // cout<< currMotifName <<"  "<< kmer[kpos] <<"  "<< currMotif(findRow(kmer[kpos]),pos+kpos) <<"\n";
        if (scoreCurrPos > MaxScore && scoreCurrPos > threshold)
          {
            MaxScore= scoreCurrPos;
            MotifName= currMotifName;
          }
      }
  }

// [[Rcpp::export]]
SEXP assign_k_to_PWMs(std::vector<std::string> kmers, Rcpp::List motifs, std::vector<std::string> motifNames, double threshold= -1000)
{
  size_t kmerNo= kmers.size();
  vector<string> motif(kmerNo);
  vector<double> scores(kmerNo);

  for (size_t k=0; k<kmerNo; k++)
  {
    string kmer= kmers[k]; string kmerRc= revComp(kmers[k]);
    string MotifName= "NA";
    double MaxScore= -100000000;
    for (size_t i = 0; i < motifs.size(); i++)
    {
      NumericMatrix currMotif= as<NumericMatrix>(motifs[i]);
      scoring_k_with_PWM(kmer,currMotif,MaxScore, MotifName, threshold, motifNames[i]);
      scoring_k_with_PWM(kmerRc,currMotif,MaxScore, MotifName, threshold, motifNames[i]);
      // cout<< motifNames[i] <<"  "<< MaxScore<<"\n";
    }
    motif[k]= MotifName;
    scores[k]= MaxScore;
    // if (k==1) stop("");
  }
  // return(fast_factor(wrap(out)));
  return(DataFrame::create(Named("motifCall")=motif, Named("maxScore")=scores ));
}






