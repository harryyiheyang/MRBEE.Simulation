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
