// [[Rcpp::plugins(cpp11)]]
// #include <Rcpp.h>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#include <iostream>
//#include <utility>
#include <cmath>
#include <vector>
#include <string>

// #include <stdlib.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

//encode binary
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

inline unsigned int  encode_oligo(const std::string& oligo, const int k, const int km1, int& currN_val){
  unsigned int  result=0;
  for (int i=0; i<k; i++) {
    result+= pow(4, km1-i) * encode_bp(oligo[i],currN_val);
  }
  return result;
}

//decode binary
inline const char decode_bp (const int x){
  switch(x){
  case 0:
    return 'A';
    break;
  case 1:
    return 'C';
    break;
  case 2:
    return 'G';
    break;
  default:
    return 'T';
  break;
  }
}
inline string decode_oligo(const unsigned int arrSlotNum, const int k, const int km1, std::string& kmerName,bool add_dot=false){
  int last2BitMask= 3;
  for (int i=0; i<k; i++) {
    kmerName[i]= decode_bp(arrSlotNum >> (2*(km1-i)) & last2BitMask);    // std::cout <<arrSlotNum<< "  >>" << (2*(km1-i))<<" "<< (arrSlotNum >> (2*(km1-i)))<< " <> " <<(arrSlotNum >> (2*(km1-i)) & last2BitMask)<<" "<< kmerName[i]<< " "<<kmerName<< "\n";
  }
  return kmerName;
}



// N is assigned with ACGT recycled
// collapse only work when all lig the same length
// [[Rcpp::export]]
SEXP kmerCntBit(std::vector<std::string> strings, int k=2, bool diffLen=false, bool collapse=false, bool asDf= true, bool all_possible_k=false, int pseudo=0)
{
  if (pseudo) all_possible_k=true; //return all kmers if using pseudo

  int kmerPerLig= strings[0].length() -k + 1;  //how many kmers per lig
  size_t ligNo = strings.size(); //total lig No.
  unsigned long kmerMask = 0;
      for (int i = 0; i < k; i++)
      {
        kmerMask <<= 2;
        kmerMask |= 3;
      }
  int km1=k-1; //for quick use
  int currN_val=0;
  uint64_t vectLen=pow(4, k);


  std::vector<int> kcnt(vectLen, pseudo); //vect for result, ini with pseudo

  std::vector< std::vector<int> > posCnt; //store pos infomation when collapse=F
  if (!collapse) {
    if (!all_possible_k)
      {posCnt.resize(vectLen);}
    else
      {posCnt.resize(vectLen, std::vector<int>(kmerPerLig,pseudo));} };

  unsigned long curr_k_code=0;
  for (size_t i = 0; i < ligNo; i++)
  {
    if (diffLen) kmerPerLig= strings[i].length() -km1;
    for (int j=0; j<kmerPerLig; j++)
    {
      if (j == 0)
        {
          curr_k_code= encode_oligo(strings[i].substr(j, k), k, km1, currN_val);
        } else
        {
          int curr_char_code = encode_bp(strings[i][j+km1], currN_val);
          curr_k_code = (curr_k_code << 2 | curr_char_code) & kmerMask;
        }

      if(!collapse)
        {
          if( (!all_possible_k) && posCnt[curr_k_code].empty()) {posCnt[curr_k_code].resize(kmerPerLig,pseudo);} //ini if not ini ed
          posCnt[curr_k_code][j]++;
        }else
        {
          kcnt[curr_k_code]++;  //#std::string sub_s= strings[i].substr(j, k); kcnt[encode_oligo(sub_s,k)]++;
        }
    }
  }



  //------------------
  std::string kmerName(k,110); //for use in decoding kmerName, ini with "n" ASCII 110
  if(all_possible_k)
    {
      vector< string >allNames(vectLen);
      for (size_t i=0; i<vectLen; i++) { allNames[i]= decode_oligo(i, k, km1,kmerName); }

      if(!collapse){
        List result= wrap(posCnt);
        result.attr("names")= allNames;
          // wrapping to df is still a bit time consuming
          IntegerVector row_name(kmerPerLig);
          for (int i=0; i< kmerPerLig; i++) {row_name[i]= i+1;} //std::to_string(i+1)
          result.attr("class") = "data.frame";
          result.attr("row.names") = row_name;
        return result;
      }else{
        if (asDf) {return DataFrame::create(Named("kmer")= allNames, Named("counts")= wrap(kcnt), _["stringsAsFactors"] = false);}
        else {IntegerVector result= wrap(kcnt); result.attr("names")= allNames; return result; }
      }
    }
        //else
        if(!collapse)
          {
            std::vector< std::vector<int> >nonEmpty;
            std::vector< std::string >nonEmpty_name;
            for (size_t i=0; i<vectLen; i++)
            {
              if (!posCnt[i].empty()) {nonEmpty.push_back(posCnt[i]); nonEmpty_name.push_back(decode_oligo(i, k, km1,kmerName)); }
            }
            List result= wrap(nonEmpty);
            result.attr("names")= nonEmpty_name;

          // wrapping to df is still a bit time consuming
            IntegerVector row_name(kmerPerLig);
            for (int i=0; i< kmerPerLig; i++) {row_name[i]= i+1;} //std::to_string(i+1)
            result.attr("class") = "data.frame";
            result.attr("row.names") = row_name;

            return result;
          }
        else  //if collapse, then consider asDF
          {
            std::vector< int >nonEmpty;
            std::vector< std::string >nonEmpty_name;
            for (size_t i=0; i<vectLen; i++)
              {
              if (kcnt[i]!=0) {nonEmpty.push_back(kcnt[i]); nonEmpty_name.push_back(decode_oligo(i, k, km1,kmerName));}
              }
            if (asDf) {return DataFrame::create(Named("kmer")= nonEmpty_name, Named("counts")= nonEmpty, _["stringsAsFactors"] = false);}
            else {IntegerVector result= wrap(nonEmpty); result.attr("names")= nonEmpty_name; return result; }
          }

}












inline const vector<unsigned long> count_k_curr (std::string& currLig, size_t &kmerPerLig,int &k,int &km1,int &currN_val, unsigned long &kmerMask)
  {
    vector<unsigned long> k_code_currLig (kmerPerLig);
    unsigned curr_k_code;
    for (size_t j=0; j<kmerPerLig; j++)
    {
      if (j == 0)
      {
        curr_k_code= encode_oligo(currLig.substr(j, k), k, km1, currN_val);
      }
      else
      {
        int curr_char_code = encode_bp(currLig[j+km1], currN_val);
        curr_k_code = (curr_k_code << 2 | curr_char_code) & kmerMask;
      }
      k_code_currLig[j]= curr_k_code;
    }
    return k_code_currLig;
  }

inline void count_gk_curr(int &k,size_t &kmerPerLig,unsigned &currLigLen,const vector<unsigned long> &k_code_currLig,vector<IntegerVector> &allkShiftCombis,IntegerVector allGapCombiLens,arma::ucube &kcnt, bool posInfo=false)
  {
    for (size_t currPos=0; currPos< kmerPerLig;currPos++) //cnt gk for each pos
    {
      size_t combiNo= allkShiftCombis.size(); //gap combi num
      size_t shiftNo= allkShiftCombis[0].size(); //gaps per combi


      for (size_t j=0; j<combiNo; j++) //loop through gap combis
      {
        unsigned long long combi_k_code=k_code_currLig[currPos];
        for (size_t idx = 0; idx < shiftNo; idx++) //pick each k segment code
        {
          combi_k_code <<= 2*k;
          if ((allGapCombiLens[j]+currPos) > currLigLen) goto end_curr_gap_combi; //if out of boundary
          combi_k_code |= k_code_currLig[allkShiftCombis[j][idx]+currPos];
        }
        kcnt(combi_k_code,j,0)++;
        if(posInfo) kcnt(combi_k_code,j,currPos+1)++; //slice[0] is for sum of all Pos
        end_curr_gap_combi:;
      }
    }
  }

bool any_sug(LogicalVector x){
  // Note the use of is_true to return a bool type
  return is_true(any(x == TRUE));
}


// N is assigned with ACGT recycled
// collapse only work when all lig the same length
// [[Rcpp::export]]
SEXP gkmerCntBit(std::vector<std::string> strings,
                 int gapNo=3,
                 int k=2,
                 IntegerVector gapMins = IntegerVector::create(2,3,0),
                 IntegerVector gapMaxs = IntegerVector::create(3,4,0),
                 int pseudo=0, bool diffLen=false, bool posInfo=false,
                 bool all_possible_k=true
)
{
  if (gapMaxs.size()!=gapNo || gapMins.size()!=gapNo) stop("gapMaxs and gapMins should be the same lens as gapNo");
  if (posInfo && diffLen) stop("cannot count with pos info when ligs has different Lengths");

  int kLenSum= (gapNo+1) * k;
  IntegerVector gapRng= (gapMaxs-gapMins+1);

  vector<IntegerVector> allGapCombis; //gen all combination of input gap
  for (int i=0; i<gapNo; i++ )
  {
    vector<int> currGaps (gapRng[i]);
    int currGapMin= gapMins[i];
    for (size_t i=0; i<currGaps.size(); i++) {currGaps[i]=currGapMin; currGapMin++;}

    vector<IntegerVector> LastCycCopy= allGapCombis;
    allGapCombis= vector< IntegerVector >();
    if (i!=0) //gen all combi
    {
      for (size_t k=0; k<LastCycCopy.size();k++)
        for (size_t j=0; j<currGaps.size();j++)
        {
          IntegerVector tmp = LastCycCopy[k];
          tmp.push_back(currGaps[j]);
          allGapCombis.push_back(tmp);
        }
    }
    else //combi of 1st gap
    {
      for (size_t i=0; i<currGaps.size(); i++) allGapCombis.push_back(IntegerVector::create(currGaps[i]));
    }
  }

  IntegerVector allGapCombiLens; //calc lens of all gap combis, stop cnt if kmer is out of lig bound
  for(size_t i=0; i<allGapCombis.size(); i++)
    allGapCombiLens.push_back(kLenSum+ sum(allGapCombis[i]));

  vector<IntegerVector> allkShiftCombis; //gen all combination of kmer shifts
  for(size_t i=0; i<allGapCombis.size(); i++)
    allkShiftCombis.push_back(cumsum(allGapCombis[i] + k));

  // counting kmers ----------------------------------------------------------------------------------------
  if (pseudo) all_possible_k=true; //return all kmers if using pseudo

  size_t kmerPerLig= strings[0].length() -k + 1;  //how many kmers per lig
  unsigned ligLens= strings[0].length();
  size_t ligNo = strings.size(); //total lig No.
  unsigned long kmerMask = 0;
  for (int i = 0; i < k; i++)
    {
      kmerMask <<= 2;
      kmerMask |= 3;
    }
  int km1=k-1; //for quick use
  int currN_val=0;

  unsigned long gapCombi= allGapCombis.size();
  unsigned long dim_k=pow(4,kLenSum);
  arma::ucube kcnt; //---------------------
    if (posInfo){ kcnt=arma::ucube(dim_k, gapCombi,ligLens+1);}
    else{kcnt=arma::ucube(dim_k, gapCombi,1);}
    kcnt.fill(pseudo);

  //counting
  for (unsigned i = 0; i < ligNo; i++)
  {
    if (diffLen) {ligLens= strings[i].length();  kmerPerLig= ligLens -km1;}
    const vector<unsigned long> k_code_currLig= count_k_curr(strings[i],kmerPerLig,k,km1,currN_val,kmerMask);
    count_gk_curr(k,kmerPerLig,ligLens,k_code_currLig,allkShiftCombis,allGapCombiLens,kcnt,posInfo);
  }



  std::vector< std::string >colNames_g(gapCombi);
  for (size_t i=0; i<gapCombi;i++)  //gen colNames
  {
    string currColName;
    for (int j =0; j< gapNo; j++) currColName+= to_string(allGapCombis[i][j])+"n";
    colNames_g[i]= currColName;
  }

  unsigned kNum= kcnt.n_rows;
  std::string kmerName(kLenSum,110); //for use in decoding kmerName, ini with "n" ASCII 110
  std::vector< std::string >rowNames_k(kNum);
  km1=kLenSum-1;

  if (!posInfo)
    {
      arma::umat kcnt_mat= kcnt.slice(0);
      if (all_possible_k) //retain all rows if all possible k
      {
        for (size_t i=0; i<kNum; i++)
          rowNames_k[i]=decode_oligo(i, kLenSum, km1,kmerName);
      }
      else
      {
        rowNames_k=std::vector< std::string >();
        uvec nonEmptyRows(kNum); unsigned nonEmptyCnt=0;
        for(size_t i=0; i<kcnt.n_rows;i++)
          if(any(kcnt_mat.row(i)!=0))
          {rowNames_k.push_back(decode_oligo(i, kLenSum, km1,kmerName)); nonEmptyRows[nonEmptyCnt]=i;nonEmptyCnt++;}
          nonEmptyRows=nonEmptyRows.head(nonEmptyCnt);
          kcnt_mat=kcnt_mat.rows(nonEmptyRows);
      }
      NumericMatrix result= wrap(kcnt_mat);
      result.attr("dimnames")= List::create(rowNames_k,colNames_g);//,seq(0,ligLens)
      return(result);
    }
  else
    {
      for (size_t i=0; i<kNum; i++)
        rowNames_k[i]=decode_oligo(i, kLenSum, km1,kmerName);

      vector<string> sliceNames;
        sliceNames.push_back("sum");
        for (unsigned i=1; i<ligLens+1;i++)
          sliceNames.push_back(to_string(i));

      List allNames=List::create(_["row"]=rowNames_k,_["col"]=colNames_g,_["slice"]=sliceNames);
      return(List::create(_["result"]=wrap(kcnt),_["dimnames"]=allNames));
    }
}



inline void scoring_gk_curr(int &k,size_t &kmerPerLig,unsigned &currLigLen,const vector<unsigned long> &k_code_currLig,vector<IntegerVector> &allkShiftCombis,IntegerVector allGapCombiLens, vector<double> &scores, arma::cube &scoreCube, bool posInfo=true)
{
  double scoreCurrLig=0;
  for (size_t currPos=0; currPos< kmerPerLig;currPos++) //cnt gk for each pos
  {
    size_t combiNo= allkShiftCombis.size(); //gap combi num
    size_t shiftNo= allkShiftCombis[0].size(); //gaps per combi


    for (size_t j=0; j<combiNo; j++) //loop through gap combis
    {
      unsigned long long combi_k_code=k_code_currLig[currPos];
      for (size_t idx = 0; idx < shiftNo; idx++) //pick each k segment code
      {
        combi_k_code <<= 2*k;
        if ((allGapCombiLens[j]+currPos) > currLigLen) goto end_curr_gap_combi; //if out of boundary
        combi_k_code |= k_code_currLig[allkShiftCombis[j][idx]+currPos];
      }

      if(posInfo) {scoreCurrLig+= scoreCube(combi_k_code,j,currPos+1);} //slice[0] is for sum of all Pos
      else{scoreCurrLig+= scoreCube(combi_k_code,j,0);}
      end_curr_gap_combi:;
    }
  }
  scores.push_back(scoreCurrLig);
}


// N is assigned with ACGT recycled
// [[Rcpp::export]]
SEXP scoring(std::vector<std::string> strings, arma::cube scoreCube,
             int gapNo=3,
             int k=2,
             IntegerVector gapMins = IntegerVector::create(2,3,0),
             IntegerVector gapMaxs = IntegerVector::create(3,4,0),
             int pseudo=0, bool diffLen=false, bool posInfo=true
            )
{
  if (gapMaxs.size()!=gapNo || gapMins.size()!=gapNo) stop("gapMaxs and gapMins should be the same lens as gapNo");
  if (posInfo && diffLen) stop("cannot count with pos info when ligs has different Lengths");

  int kLenSum= (gapNo+1) * k;
  IntegerVector gapRng= (gapMaxs-gapMins+1);

  vector<IntegerVector> allGapCombis; //gen all combination of input gap
  for (int i=0; i<gapNo; i++ )
  {
    vector<int> currGaps (gapRng[i]);
    int currGapMin= gapMins[i];
    for (size_t i=0; i<currGaps.size(); i++) {currGaps[i]=currGapMin; currGapMin++;}

    vector<IntegerVector> LastCycCopy= allGapCombis;
    allGapCombis= vector< IntegerVector >();
    if (i!=0) //gen all combi
    {
      for (size_t k=0; k<LastCycCopy.size();k++)
        for (size_t j=0; j<currGaps.size();j++)
        {
          IntegerVector tmp = LastCycCopy[k];
          tmp.push_back(currGaps[j]);
          allGapCombis.push_back(tmp);
        }
    }
    else //combi of 1st gap
    {
      for (size_t i=0; i<currGaps.size(); i++) allGapCombis.push_back(IntegerVector::create(currGaps[i]));
    }
  }

  IntegerVector allGapCombiLens; //calc lens of all gap combis, stop cnt if kmer is out of lig bound
  for(size_t i=0; i<allGapCombis.size(); i++)
    allGapCombiLens.push_back(kLenSum+ sum(allGapCombis[i]));

  vector<IntegerVector> allkShiftCombis; //gen all combination of kmer shifts
  for(size_t i=0; i<allGapCombis.size(); i++)
    allkShiftCombis.push_back(cumsum(allGapCombis[i] + k));


  // counting kmers ----------------------------------------------------------------------------------------
  // if (pseudo) all_possible_k=true; //return all kmers if using pseudo

  size_t kmerPerLig= strings[0].length() -k + 1;  //how many kmers per lig
  unsigned ligLens= strings[0].length();         //cout << ligLens << "  " << scoreCube.n_slices <<"\n";
    if ((scoreCube.n_cols!=allGapCombis.size())) stop("gaps of target lib different from that of the scoreCube");
      else if (scoreCube.n_rows!=pow(4, kLenSum)) stop("kmerLens of target lib different from that of the scoreCube");
      else if (scoreCube.n_slices!=(ligLens+1)) stop("slices of target lib different from that of the scoreCube");

  size_t ligNo = strings.size(); //total lig No.
  unsigned long kmerMask = 0;
    for (int i = 0; i < k; i++)
    {
      kmerMask <<= 2;
      kmerMask |= 3;
    }
  int km1=k-1; //for quick use
  int currN_val=0;

  vector<double> scores; //---------------------

  //counting
  for (unsigned i = 0; i < ligNo; i++)
  {
    // if (diffLen) {ligLens= strings[i].length();  kmerPerLig= ligLens -km1;}
    const vector<unsigned long> k_code_currLig= count_k_curr(strings[i],kmerPerLig,k,km1,currN_val,kmerMask);
    scoring_gk_curr(k,kmerPerLig,ligLens,k_code_currLig,allkShiftCombis,allGapCombiLens,scores,scoreCube,posInfo);
  }
  return wrap(scores);
}







//#include  <RInside.h>                   // for the embedded R via RInside
//
//int main(int argc, char *argv[]) {
//
//    // create an embedded R instance
//    RInside R(argc, argv);
//
//    // demonstrates setting R variables to C++ variables
//    std::string outPath = "~/Desktop";
//    std::string outFile = "myPlot.png";
//    R["outPath"] = outPath;
//    R["outFile"] = outFile;
//
//    // build the sequence of R commands as a string
//
//    std::string cmd ="library(readr); seq=read_csv('~/Nut_zhuData/Analysis/Analysis2/kcnt_test/testseq.txt',col_names = F); seq$X1 ";
//    std::vector<std::string> seq= R.parseEval(cmd);
//
//    SEXP kcnt2= kmerCntBit(seq, 2, false, true,true, false, 2);
//    return 0;
//}


/*** R

*/
