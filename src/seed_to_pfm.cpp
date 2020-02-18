// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <iostream>
// #include "common_op.h"
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
NumericMatrix matRevComp(NumericMatrix seqs)
{
  NumericMatrix rc(seqs.nrow(),seqs.ncol());  //NumericMatrix(seqs.nrow(),seqs.ncol());
  for (size_t i=0; i<seqs.size(); i++) rc[i]= seqs[i];
  std::reverse(rc.begin(),rc.end());
  // std::reverse(seqs.begin(),seqs.end());
  return rc;
}


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

std::string revComp1(std::string seq)
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

// only allow hamming dist 1
// [[Rcpp::export]]
NumericMatrix pfm_from_seed_notUsed(std::vector<std::string> seqs, std::string seed1="AAA",int gapLen=1, std::string seed2="AAA",
                            NumericVector seed1_start= NumericVector::create(2,3,4), int flankLen=4, bool all_start_with_specified_gap=true) {
  int seed1_len=seed1.length(); int seed2_len= seed2.length();
  int seq_length=seqs[0].size();
  NumericMatrix pfm(4,seed1_len+gapLen+seed2_len+(2*flankLen));
  std::fill(pfm.begin(),pfm.end(),1);
  seed1_start= seed1_start-1; //corr to start from 0
  if (all_start_with_specified_gap) // add all possible start positions
  {
    // cout<<"all possible start pos!!\n";
    seed1_start=NumericVector::create();
    for (size_t i=0; i<(seq_length-seed1_len-gapLen-seed2_len+1); i++) seed1_start.push_back(i);
  }
  NumericVector seed2_start= seed1_start + seed1_len + gapLen;
  NumericVector flank1_start=seed1_start-flankLen;
  NumericVector flank2_start=seed2_start+seed2_len;
  NumericVector gap_start=seed1_start+seed1_len;
  int currN_val=0; int seed1_pfmstart=flankLen; int seed2_pfmstart=flankLen+seed1_len+gapLen;
  int gap_pfmstart=flankLen+seed1_len; int flank2_pfmstart=flankLen+seed1_len+gapLen+seed2_len;

  // , bool two_strands=false /// add in param/////
  // if (two_strands)
  // {
  //   // cout<<"2 strand!!\n";
  //   seqs.resize(seqs.size()*2);
  //   int half_size=seqs.size()/2;
  //   for (size_t i=half_size; i<seqs.size(); i++) seqs[i]=revComp1(seqs[i-half_size]);
  // }


  for (size_t i=0; i<seqs.size(); i++)
  {
    for (size_t j=0; j<seed1_start.size();j++) //loop through all start pos
    {
      std::string str1= seqs[i].substr(seed1_start[j],seed1_len);
      std::string str2= seqs[i].substr(seed2_start[j],seed2_len);

      int HammingDist=0;
      int seed1_diff_pos=-1; int seed2_diff_pos=-1;
      //dist of str1
      for (size_t ind=0; ind<seed1_len; ind++)
      {
        if (str1[ind] != seed1[ind])
        {
          if (HammingDist>0) goto outer;
          HammingDist++; seed1_diff_pos=ind;
        }
      }
      //dist of str2
      for (size_t ind=0; ind<seed2_len; ind++)
      {
        if (str2[ind] != seed2[ind])
        {
          if (HammingDist>0) goto outer;
          HammingDist++; seed2_diff_pos=ind;
        }
      }


      // start counting if dist <=1
      if (HammingDist==0) //count flank
      {
        //count seed pos
        for (size_t ind=0; ind<seed1_len; ind++)
          pfm(encode_bp(str1[ind],currN_val),seed1_pfmstart+ind)++;
        for (size_t ind=0; ind<seed2_len; ind++)
          pfm(encode_bp(str2[ind],currN_val),seed2_pfmstart+ind)++;

        //count flanking pos, only if flanks is not out of bound
        if(flank1_start[j]>0)
        {
          std::string flank1= seqs[i].substr(flank1_start[j],flankLen);
          for (size_t ind=0; ind<flankLen; ind++)
            pfm(encode_bp(flank1[ind],currN_val),ind)++;
        }
        if(flank2_start[j]+flankLen < (seq_length-1))
        {
          std::string flank2= seqs[i].substr(flank2_start[j],flankLen);
          for (size_t ind=0; ind<flankLen; ind++)
            pfm(encode_bp(flank2[ind],currN_val),flank2_pfmstart+ind)++;
        }

        //count gap pos
        std::string gap= seqs[i].substr(gap_start[j],gapLen);
        for (size_t ind=0; ind<gapLen; ind++)
          pfm(encode_bp(gap[ind],currN_val), gap_pfmstart+ind)++;
      }
      else if (HammingDist==1) //count seed pos only
      {
        if (seed1_diff_pos!=-1){ pfm(encode_bp(str1[seed1_diff_pos],currN_val),seed1_pfmstart+seed1_diff_pos)++; }
        else if(seed2_diff_pos!=-1){ pfm(encode_bp(str2[seed2_diff_pos],currN_val),seed2_pfmstart+seed2_diff_pos)++; }
        else{stop("seed1_diff_pos or seed2_diff_pos should not all be -1");}
      }
      else {stop("Hamming Dist is not supposted to be >1");}

      outer:;
    }

  }
  return(pfm);
}



// // only allow hamming dist 1
// // [[Rcpp::export]]
// NumericMatrix pfm_from_seed_IUPAC_notUsed(std::vector<std::string> seqs, std::string seed="AAA",
//                                     NumericVector seed_start= NumericVector::create(2,3,4), bool all_start_with_specified_gap=true) {
//   int seed_len=seed.length();
//   int seq_length=seqs[0].size();
//   NumericMatrix pfm(4,seed_len);
//   std::fill(pfm.begin(),pfm.end(),1);
//   seed_start= seed_start-1; //corr to start from 0
//   if (all_start_with_specified_gap) // add all possible start positions
//   {
//     // cout<<"all possible start pos!!\n";
//     seed_start=NumericVector::create();
//     for (size_t i=0; i<(seq_length-seed_len+1); i++) seed_start.push_back(i);
//   }
//
//   int currN_val=0; int seed1_pfmstart=flankLen; int seed2_pfmstart=flankLen+seed1_len+gapLen;
//   int gap_pfmstart=flankLen+seed1_len; int flank2_pfmstart=flankLen+seed1_len+gapLen+seed2_len;
//
//   // , bool two_strands=false /// add in param/////
//   // if (two_strands)
//   // {
//   //   // cout<<"2 strand!!\n";
//   //   seqs.resize(seqs.size()*2);
//   //   int half_size=seqs.size()/2;
//   //   for (size_t i=half_size; i<seqs.size(); i++) seqs[i]=revComp1(seqs[i-half_size]);
//   // }
//
//
//   for (size_t i=0; i<seqs.size(); i++)
//   {
//     for (size_t j=0; j<seed1_start.size();j++) //loop through all start pos
//     {
//       std::string str1= seqs[i].substr(seed1_start[j],seed1_len);
//       std::string str2= seqs[i].substr(seed2_start[j],seed2_len);
//
//       int HammingDist=0;
//       int seed1_diff_pos=-1; int seed2_diff_pos=-1;
//       //dist of str1
//       for (size_t ind=0; ind<seed1_len; ind++)
//       {
//         if (str1[ind] != seed1[ind])
//         {
//           if (HammingDist>0) goto outer;
//           HammingDist++; seed1_diff_pos=ind;
//         }
//       }
//       //dist of str2
//       for (size_t ind=0; ind<seed2_len; ind++)
//       {
//         if (str2[ind] != seed2[ind])
//         {
//           if (HammingDist>0) goto outer;
//           HammingDist++; seed2_diff_pos=ind;
//         }
//       }
//
//
//       // start counting if dist <=1
//       if (HammingDist==0) //count flank
//       {
//         //count seed pos
//         for (size_t ind=0; ind<seed1_len; ind++)
//           pfm(encode_bp(str1[ind],currN_val),seed1_pfmstart+ind)++;
//         for (size_t ind=0; ind<seed2_len; ind++)
//           pfm(encode_bp(str2[ind],currN_val),seed2_pfmstart+ind)++;
//
//         //count flanking pos, only if flanks is not out of bound
//         if(flank1_start[j]>0)
//         {
//           std::string flank1= seqs[i].substr(flank1_start[j],flankLen);
//           for (size_t ind=0; ind<flankLen; ind++)
//             pfm(encode_bp(flank1[ind],currN_val),ind)++;
//         }
//         if(flank2_start[j]+flankLen < (seq_length-1))
//         {
//           std::string flank2= seqs[i].substr(flank2_start[j],flankLen);
//           for (size_t ind=0; ind<flankLen; ind++)
//             pfm(encode_bp(flank2[ind],currN_val),flank2_pfmstart+ind)++;
//         }
//
//         //count gap pos
//         std::string gap= seqs[i].substr(gap_start[j],gapLen);
//         for (size_t ind=0; ind<gapLen; ind++)
//           pfm(encode_bp(gap[ind],currN_val), gap_pfmstart+ind)++;
//       }
//       else if (HammingDist==1) //count seed pos only
//       {
//         if (seed1_diff_pos!=-1){ pfm(encode_bp(str1[seed1_diff_pos],currN_val),seed1_pfmstart+seed1_diff_pos)++; }
//         else if(seed2_diff_pos!=-1){ pfm(encode_bp(str2[seed2_diff_pos],currN_val),seed2_pfmstart+seed2_diff_pos)++; }
//         else{stop("seed1_diff_pos or seed2_diff_pos should not all be -1");}
//       }
//       else {stop("Hamming Dist is not supposted to be >1");}
//
//       outer:;
//     }
//
//   }
//   return(pfm);
// }



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/

