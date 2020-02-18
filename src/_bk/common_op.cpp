// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <iostream>
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


inline const short int encode_bp (const char x, int& currN_val){
  switch(x){
  case 'N':case 'n': //illegal value
    return (currN_val++ % 4);
  case 'A':case 'a':
    return 0;
    break;
  case 'C':case 'c':
    return 1;
    break;
  case 'G':case 'g':
    return 2;
    break;
  default:
    return 3;
  break;
  }
}
