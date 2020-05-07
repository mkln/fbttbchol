#include "RcppArmadillo.h"

arma::mat ec(const arma::mat& x, const arma::mat& y, double sigmasq, double phi, bool same=false){
  // 0 based indexing
  if(same){
    arma::mat pmag = arma::sum(x % x, 1);
    int np = x.n_rows;
    arma::mat K = sigmasq * exp(-phi * sqrt(abs(arma::repmat(pmag.t(), np, 1) + arma::repmat(pmag, 1, np) - 2 * x * x.t())));
    return K;
  } else {
    arma::mat pmag = arma::sum(x % x, 1);
    arma::mat qmag = arma::sum(y % y, 1);
    int np = x.n_rows;
    int nq = y.n_rows;
    arma::mat K = sigmasq * exp(-phi * sqrt(abs(arma::repmat(qmag.t(), np, 1) + arma::repmat(pmag, 1, nq) - 2 * x * y.t())));
    return K;
  }
}

//[[Rcpp::export]]
arma::field<arma::mat> make_covlist(const arma::mat& coords, int nsim=10){
  arma::field<arma::mat> C(nsim);
  
#pragma omp parallel for
  for(int i=0; i<nsim; i++){
    double phi = arma::conv_to<double>::from( exp( arma::randn(1) ));
    C(i) = ec(coords, coords, 1, phi, true);
  }
  return C;
}


extern "C" {
  void mb02cdmod_( char* job, char* typet, int* K, int* N, double* T, int* LDT, 
                   double* G, int* LDG, double* R, int* LDR, double* L, int* LDL,
                   double* CS, int* LCS, double* DWORK, int* LDWORK, int* INFO );
}


Rcpp::List r_mb02cd(const arma::mat& C, int k, int n, char jobin) {
  arma::mat T = C.cols(0, k-1);
  
  char job(jobin);
  char typet('C');
  
  int ldt = n*k;
  int ldg = n*k;
  int ldr = n*k;
  int ldl = n*k;
  int lcs = 3*(n-1)*k;
  int ldwork = 3*n*k;
  
  arma::mat G = arma::zeros(ldg, 2*k);
  arma::mat R = arma::zeros(ldr, n*k);
  arma::mat L = arma::zeros(ldl, n*k);
  arma::vec CS = arma::zeros(lcs);
  arma::vec DWORK = arma::zeros(ldwork);
  
  int info=99;
  
  mb02cdmod_(&job, &typet, &k, &n, 
          T.memptr(), &ldt,
          G.memptr(), &ldg,
          R.memptr(), &ldr,
          L.memptr(), &ldl,
          CS.memptr(), &lcs,
          DWORK.memptr(), &ldwork,
          &info);
  
  Rcpp::Rcout << info << "\n";
  
  return Rcpp::List::create(
    Rcpp::Named("T") = T,
    Rcpp::Named("G") = G,
    Rcpp::Named("R") = R,
    Rcpp::Named("L") = L,
    Rcpp::Named("CS") = CS,
    Rcpp::Named("DWORK") = DWORK,
    Rcpp::Named("info") = info
  );
}

//[[Rcpp::export]]
arma::mat arma_chol(const arma::mat& x){
  return arma::chol(x, "lower");
}

//[[Rcpp::export]]
arma::mat arma_invchol(const arma::mat& x){
  return arma::inv(arma::trimatl(arma::chol(x, "lower")));
}

//[[Rcpp::export]]
arma::field<arma::mat> arma_invchol_list(const arma::mat& coords, const arma::field<arma::mat>& covlist,
                                         bool print = false){
  int nsim = covlist.n_elem;
  arma::field<arma::mat> result(nsim);
  
  std::chrono::steady_clock::time_point start, end;
  start = std::chrono::steady_clock::now();
  
#pragma omp parallel for
  for(int i=0; i<nsim; i++){
    result(i) = arma::inv(arma::trimatl(arma::chol(covlist(i), "lower")));
  }
  end = std::chrono::steady_clock::now();
  double timing = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  if(print){
    Rcpp::Rcout << timing << "ms\n";
  }
  
  return result;
}



//[[Rcpp::export]]
arma::mat f_chol(const arma::mat& C, int k, int n) {
  // chol of BTTB
  char job('R');
  char typet('C');
  
  int ldt = n*k;
  int ldg = n*k;
  int ldr = n*k;
  int ldl = n*k;
  int lcs = 3*(n-1)*k;
  int ldwork = n*k;
  
  arma::mat G = arma::zeros(ldg, 2*k);
  arma::mat R = arma::zeros(ldr, n*k);
  arma::mat L = arma::zeros(ldl, n*k);
  arma::vec CS = arma::zeros(lcs);
  arma::vec DWORK = arma::zeros(ldwork);
  int info=99;
  
  arma::mat T = C.cols(0, k-1);
  mb02cdmod_(&job, &typet, &k, &n, 
             T.memptr(), &ldt,
             G.memptr(), &ldg,
             R.memptr(), &ldr,
             L.memptr(), &ldl,
             CS.memptr(), &lcs,
             DWORK.memptr(), &ldwork,
             &info);
  return R;
}

//[[Rcpp::export]]
arma::mat f_invchol(const arma::mat& C, int k, int n) {
  // inverse chol of BTTB
  char job('L');
  char typet('C');
  
  int ldt = n*k;
  int ldg = n*k;
  int ldr = n*k;
  int ldl = n*k;
  int lcs = 3*(n-1)*k;
  int ldwork = n*k;
  
  arma::mat G = arma::zeros(ldg, 2*k);
  arma::mat R = arma::zeros(ldr, n*k);
  arma::mat L = arma::zeros(ldl, n*k);
  arma::vec CS = arma::zeros(lcs);
  arma::vec DWORK = arma::zeros(ldwork);
  int info=99;
  
  arma::mat T = C.cols(0, k-1);
  mb02cdmod_(&job, &typet, &k, &n, 
             T.memptr(), &ldt,
             G.memptr(), &ldg,
             R.memptr(), &ldr,
             L.memptr(), &ldl,
             CS.memptr(), &lcs,
             DWORK.memptr(), &ldwork,
             &info);
  return L.t();
}


//[[Rcpp::export]]
arma::field<arma::mat> f_invchol_list(const arma::mat& coords, int k, int n, const arma::field<arma::mat>& covlist,
                                      bool print=false) {
  // inverse of chol of BTTB (on a list of matrices)
  char job('L');
  char typet('C');
  
  int nsim = covlist.n_elem;
  
  arma::field<arma::mat> result(nsim);
  
  std::chrono::steady_clock::time_point start, end;
  start = std::chrono::steady_clock::now();
  
#pragma omp parallel for
  for(int i=0; i<nsim; i++){
    arma::mat T = covlist(i).cols(0, k-1);
  
    int ldt = n*k;
    int ldg = n*k;
    int ldr = n*k;
    int ldl = n*k;
    int lcs = 3*(n-1)*k;
    int ldwork = n*k;
    
    arma::mat G = arma::zeros(ldg, 2*k);
    arma::mat R = arma::zeros(ldr, n*k);
    arma::mat L = arma::zeros(ldl, n*k);
    arma::vec CS = arma::zeros(lcs);
    arma::vec DWORK = arma::zeros(ldwork);
    int info=99;

    mb02cdmod_(&job, &typet, &k, &n, 
            T.memptr(), &ldt,
            G.memptr(), &ldg,
            R.memptr(), &ldr,
            L.memptr(), &ldl,
            CS.memptr(), &lcs,
            DWORK.memptr(), &ldwork,
            &info);
    result(i) = L.t();
  }
  end = std::chrono::steady_clock::now();
  double timing = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  if(print){
    Rcpp::Rcout << timing << "ms\n";
  }
  return result;
}
