//counting kmers in R
//use:
//kmerCnt(strArr, k=2, collapse=F)
//NB!!!!!!!!!!!!!!!  each kmer should be the same length !!!!!!!!!!!!!!!!!!!!!!!!!


// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <map>
#include <string.h>
#include <unordered_map>
// #include <iostream>

using namespace Rcpp;
// using namespace std;

// [[Rcpp::export]]
SEXP kmerCnt(std::vector< std::string > strings, int k=2, bool collapse=false, bool diffLen=false, bool asDf= false)
{
  int kmerNo= strings[0].length() -k + 1;  //how many kmers per lig
  int ligNo = strings.size();

  if (!collapse)  //retain pos info if not collapse
  {
      // storage of the out put
      std::map<std::string, std::vector<int> > kcnt;

      for (int i = 0; i < ligNo; i++) // each lig
      {
        if (diffLen) //decide params for each lig if they differ in length
        {
          kmerNo= strings[i].length() -k + 1;
        }
        for (int j=0; j<kmerNo; j++) // each pos of the win
        {
          std::string sub_s= strings[i].substr(j, k);
          if (  strchr(sub_s.c_str(), 'N') == NULL && strchr(sub_s.c_str(), 'n') == NULL )
          {
            if (kcnt[sub_s].empty()) {kcnt[sub_s].resize(kmerNo);}
            kcnt[sub_s][j]++;
          }
        }
      }
      return Rcpp::DataFrame(kcnt);
  }
  else  //collapse if not require pos info
  {
    // storage of the out put
    std::map<std::string, int> kcnt;

    for (int i = 0; i < ligNo; i++) // each lig
    {
      if (diffLen) //decide params for each lig if they differ in length
      {
        kmerNo= strings[i].length() -k + 1;
      }
      for (int j=0; j<kmerNo; j++) // each pos of the win
      {
        std::string sub_s= strings[i].substr(j, k);
        if (  strchr(sub_s.c_str(), 'N') == NULL && strchr(sub_s.c_str(), 'n') == NULL )
        {
          // if (kcnt[sub_s].empty()) {kcnt[sub_s].resize(kmerNo);}
          kcnt[sub_s]++;
        }
      }
    }

    if (asDf) //convert result to df
    {
      int keyNo= kcnt.size();
      CharacterVector kfull(keyNo);
      IntegerVector counts(keyNo);

      typedef std::map<std::string, int>::iterator it_type; int cnt=0;
      for(it_type iterator = kcnt.begin(); iterator != kcnt.end(); iterator++)
      {
        kfull[cnt]=iterator->first;
        counts[cnt]=iterator->second;
        cnt++;
      }
      return DataFrame::create(Named("kmer")= kfull, Named("counts")= counts, _["stringsAsFactors"] = false);
    }

    return wrap(kcnt);
  }

}


// // [[Rcpp::export]]
// SEXP kmerCnt_gap(std::vector< std::string > strings, int k=2, int gapLen=3)
// {
//   int kmerNo= strings[0].length() -k*2 - gapLen + 1;  //how many kmers per lig
//   int ligNo = strings.size();
//   std::string gapLen_str= "X"+std::to_string(gapLen)+"X";
//
//     // storage of the out put
//     std::map<std::string, int> kcnt;
//
//     for (int i = 0; i < ligNo; i++) // each lig
//     {
//       for (int j=0; j<kmerNo; j++) // each pos of the win
//       {
//         std::string sub_s= strings[i].substr(j, k) + gapLen_str+ strings[i].substr(j+k+gapLen, k);
//         if (  strchr(sub_s.c_str(), 'N') == NULL && strchr(sub_s.c_str(), 'n') == NULL )
//         {
//           // if (kcnt[sub_s].empty()) {kcnt[sub_s].resize(kmerNo);}
//           kcnt[sub_s]++;
//         }
//       }
//
//     }
//     return wrap(kcnt);
// }


// count kmers with all gapLen
// [[Rcpp::export]]
SEXP kmerCnt_allgap(std::vector< std::string > strings, int k=2, int maxGap=0, bool asDf=true, bool diffLen=false, int minGap=0) // min/maxGap=0 if no limit
{
  int ligLen=strings[0].length();
  int kmerNo_nogap= ligLen -k*2 + 1;  //how many kmers per lig
  int maxGapLen=ligLen- k*2;
  int ligNo = strings.size();

  std::vector<std::string> gapLen_str; gapLen_str.resize(maxGapLen + 1);
  for (int gapLen=0; gapLen<=maxGapLen; gapLen++)
      gapLen_str[gapLen]= "X"+std::to_string(gapLen)+"X";

  // storage of the out put
  std::map<std::string, int> kcnt;

  for (int i = 0; i < ligNo; i++) // each lig
  {
        if (diffLen) //decide params for each lig if they differ in length
          {
            ligLen=strings[i].length();
            kmerNo_nogap= ligLen -k*2 + 1;  //how many kmers per lig
            maxGapLen=ligLen- k*2;
          }
        for (int j=0, maxGapLenCurr=maxGapLen; j<kmerNo_nogap; j++, maxGapLenCurr--) // each pos of the win
        {
          std::string sub_s1= strings[i].substr(j, k);
          for (int gapLen=0; gapLen<=maxGapLenCurr; gapLen++)
          {
            if(maxGap && gapLen>maxGap){break;}
            if(minGap && gapLen<minGap){continue;}
            std::string sub_s= sub_s1 + gapLen_str[gapLen]+ strings[i].substr(j+k+gapLen, k);
            if (  strchr(sub_s.c_str(), 'N') == NULL && strchr(sub_s.c_str(), 'n') == NULL )
            {
              kcnt[sub_s]++;
            }
          }
        }
  }


  if (asDf) //convert result to df
  {
    int keyNo= kcnt.size();
    CharacterVector kfull(keyNo);
    CharacterVector k1(keyNo);
    IntegerVector gapLen(keyNo);
    CharacterVector k2(keyNo);
    IntegerVector counts(keyNo);

    typedef std::map<std::string, int>::iterator it_type; int cnt=0;
    for(it_type iterator = kcnt.begin(); iterator != kcnt.end(); iterator++)
    {
      std::string kmer= iterator->first;
      kfull[cnt]=kmer;
      k1[cnt]=kmer.substr(0,k);
      gapLen[cnt]=atoi(kmer.substr(k+1,kmer.length()-k-k-2).c_str());
      k2[cnt]=kmer.substr(kmer.length()-k,k);
      counts[cnt]=iterator->second;
      cnt++;
    }
    return DataFrame::create(_["kmer"]= kfull, _["kmer1"]= k1, _["gapLen"]= gapLen, _["kmer2"]= k2, _["counts"]= counts,_["stringsAsFactors"] = false);
  }

  return wrap(kcnt);
}





/*** R

*/
