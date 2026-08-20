#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
S4 sddmm_csc(S4 Y, NumericMatrix A, NumericMatrix B) {
    IntegerVector p = Y.slot("p");
    IntegerVector i = Y.slot("i");
    IntegerVector Dim = Y.slot("Dim");

    int M = Dim[1];
    int K = A.ncol();
    int nnz = i.size();

    NumericVector out_x(nnz);

    for (int j = 0; j < M; j++) {
        for (int idx = p[j]; idx < p[j + 1]; idx++) {
            int row = i[idx];
            double dot = 0.0;
            for (int k = 0; k < K; k++) {
                dot += A(row, k) * B(k, j);
            }
            out_x[idx] = dot;
        }
    }

    S4 result("dgCMatrix");
    result.slot("i") = i;
    result.slot("p") = p;
    result.slot("x") = out_x;
    result.slot("Dim") = Dim;
    result.slot("Dimnames") = Y.slot("Dimnames");

    return result;
}
