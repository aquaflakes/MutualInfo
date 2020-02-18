// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
// #include <RcppArmadillo.h>

#include <iostream>
// #include <vector>
// #include <array>

using namespace Rcpp;
using namespace std;


// SEXP readCube(SEXP myArray)
// {
//   NumericVector vecArray(myArray);
//   IntegerVector arrayDims = vecArray.attr("dim");
//
//   arma::cube cubeArray(vecArray.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
//
//   //change one element in the array/cube
//   cubeArray(0,0,0) = 518;
//
//   return(wrap(cubeArray));
// }


NumericMatrix row_erase (NumericMatrix& x, IntegerVector& rowID) {
  rowID = rowID.sort();

  NumericMatrix x2(Dimension(x.nrow()- rowID.size(), x.ncol()));
  int iter = 0;
  int del = 1; // to count deleted elements
  for (int i = 0; i < x.nrow(); i++) {
    if (i != rowID[del - 1]) {
      x2.row(iter) = x.row(i);
      iter++;
    } else {
      del++;
    }
  }
  return x2;
}

NumericMatrix col_erase (NumericMatrix& x, IntegerVector& colID) {
  colID = colID.sort();

  NumericMatrix x2(Dimension(x.nrow(), x.ncol()- colID.size()));
  int iter = 0;
  int del = 1;
  for (int i = 0; i < x.ncol(); i++) {
    if (i != colID[del - 1]) {
      x2.column(iter) = x.column(i);
      iter++;
    } else {
      del++;
    }
  }
  return x2;
}



IntegerVector cumsum1(IntegerVector x){
  // initialize an accumulator variable
  double acc = 0;

  // initialize the result vector
  IntegerVector res(x.size());

  for(int i = 0; i < x.size(); i++){
    acc += x[i];
    res[i] = acc;
  }
  return res;
}
