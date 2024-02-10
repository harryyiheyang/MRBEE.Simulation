require(Rcpp)
require(RcppArmadillo)

cppFunction(depends = "RcppArmadillo", code = '
            arma::mat matrixInverse(const arma::mat& A) {
            return inv(A);
            }
            ')

cppFunction(depends = "RcppArmadillo", code = '
            arma::mat matrixMultiply(const arma::mat& A, const arma::mat& B) {
            return A * B;
            }
            ')

Rcpp::cppFunction(depends = "RcppArmadillo", code = '
arma::vec matrixVectorMultiply(const arma::mat& A, const arma::vec& b) {
  return A * b;
}
')

cppFunction(depends = "RcppArmadillo", code = '
arma::mat kroneckerProduct(const arma::mat& A, const arma::mat& B) {
  arma::mat K(A.n_rows * B.n_rows, A.n_cols * B.n_cols);
  for (std::size_t i = 0; i < A.n_rows; ++i) {
    for (std::size_t j = 0; j < A.n_cols; ++j) {
      K.submat(i * B.n_rows, j * B.n_cols, (i + 1) * B.n_rows - 1, (j + 1) * B.n_cols - 1) = A(i, j) * B;
    }
  }
  return K;
}
')

cppFunction(depends = "RcppArmadillo", code = '
arma::mat Corr(const arma::mat& A) {
  int n = A.n_rows;
  
  // Remove the mean of each column
  arma::mat B = A.each_row() - mean(A, 0);
  
  // Compute the sample covariance matrix
  arma::mat S = trans(B) * B / (n - 1);
  
  // Convert the sample covariance matrix to the sample correlation matrix
  arma::vec D = 1.0 / sqrt(S.diag());
  D.replace(arma::datum::nan, 0);  // Replace NaN values with 0
  arma::mat R = S.each_col() % D;
  R = R.each_row() % trans(D);
  
  return R;
}
')
