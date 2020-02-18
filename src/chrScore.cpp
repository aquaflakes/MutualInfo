//kmer scores for each Chr position

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <iostream>
#include <cmath>

#include <vector>
#include <string>
#include <map>

// #include <set>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::NumericVector getChrScore_scored_k( std::vector<std::string> kmers, std::vector<double> scores, std::string chrSeq, int kmerLen)
{

  std::map<std::string, double> topKmers;
  for(size_t i=0; i< kmers.size(); i++ )
    {
      topKmers[kmers[i]]= scores[i];
    }

  int kmerNo = chrSeq.length()- kmerLen + 1;
  vector<double> score (kmerNo,0);


  for (int i=0; i<kmerNo; i++) // each pos of the win
  {
    std::string curr_k= chrSeq.substr(i, kmerLen);
    if (topKmers.find(curr_k) != topKmers.end()) {score[i]= topKmers[curr_k];}
  }

  return wrap(score);
}


// [[Rcpp::export]]
Rcpp::IntegerVector getChrScore(std::vector< std::string > topKmers, std::string chrSeq)
{
  int kmerLen= topKmers[0].length();  //how many kmers per lig
  int kmerNo = chrSeq.length()- kmerLen + 1;

  std::set<std::string> topKmerSet (topKmers.begin(), topKmers.end());
  vector<int> score (kmerNo,0);


  for (int i=0; i<kmerNo; i++) // each pos of the win
  {
    if (topKmerSet.find(chrSeq.substr(i, kmerLen)) != topKmerSet.end()) score[i]++;
  }

  return wrap(score);
}




/*** R

*/
