
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
CharacterMatrix seqFregments(std::vector<std::string> strings, int k=2L)
{
  size_t kmerPerLig= strings[0].length() -k + 1;
  size_t ligNo = strings.size();
  CharacterMatrix result(ligNo,kmerPerLig);

  for (size_t i = 0; i < ligNo; i++)
    for (int j=0; j<kmerPerLig; j++)
      result(i,j)= strings[i].substr(j, k);
  return wrap(result);
}
