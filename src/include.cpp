#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::mat covmat_c(arma::mat Xt, arma::vec mu, arma::vec u){
  int P = Xt.n_cols;
  int TIME = Xt.n_rows;
  arma::mat covest(P,P); covest.fill(0);
  for(int time = 0; time < TIME; time++){
    arma::vec Xtmu = Xt.row(time).t() - mu;
    covest = covest + u(time)*Xtmu*Xtmu.t();
  }
  covest = covest/accu(u);
  return(covest);
}

// [[Rcpp::export]]
arma::mat standardize_rows_cpp(arma::mat A){
  arma::mat Anew = A;
  for(int i = 0; i < A.n_rows; i++){
    arma::rowvec A_i = A.row(i);
    Anew.row(i) = A_i/accu(A_i.t());
  }
  return(Anew);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnrm_arma(const arma::mat& x,
                      const arma::rowvec& mean,
                      const arma::mat& sigma,
                      bool logd = false,
                      bool chol = false) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat root(sigma.n_rows,sigma.n_cols);
  // Check if the input was actually the inverse square root precomputed
  if(chol){
    root = sigma;
  }
  else{
    root = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  }
  double rootsum = arma::sum(log(root.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = root * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootsum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
arma::vec logvec_c(const arma::vec& x){
  arma::vec result(x.n_elem); result.fill(0);
  for(int i = 0; i < x.n_elem; i++){
    result(i) = log(x(i));
  }
  return(result);
}

// [[Rcpp::export]]
double logsumexp_cpp(arma::vec X){
  double maxvec = X.max();
  return(maxvec + log(accu(exp(X-maxvec))));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::umat get_locations(arma::sp_mat& B)
{
  // Obtains a list of coordinates from the sparse matrix by using iterators
  // First calculates total number of points and, then, obtains list.
  // Reference: https://stackoverflow.com/questions/40222092/access-and-modify-the-non-zero-elements-of-sparse-matrix-of-class-armasp-mat-u

  // Make const iterator
  arma::sp_mat::const_iterator start = B.begin();
  arma::sp_mat::const_iterator end   = B.end();

  // Calculate number of points
  int n = std::distance(start, end);

  // Kill process if no values are found (very sparse matrix)
  if (n <= 0) { Rcpp::stop("No values found!"); }

  // Build a location storage matrix
  arma::umat locs(2, n);

  // Create a vector to store each row information in. (Row, Col)
  arma::uvec temp(2);

  // Start collecting locations
  arma::sp_mat::const_iterator it = start;
  for(int i = 0; i < n; ++i)
  {
    temp(0) = it.row();
    temp(1) = it.col();
    locs.col(i) = temp;
    ++it; // increment
  }

  return locs;
}

// [[Rcpp::export]]
void vec2params_covar(const arma::vec& params_vec, int num_states, int Z_n_cols,
                      arma::mat& betahat){
  // Get the covariate effects
  int index = 0;
  for(int i = 0; i < num_states; i++){
    for(int j = 0; j < Z_n_cols; j++){
      betahat(i,j) = params_vec(index);
      index++;
    }
  }
}


// [[Rcpp::export]]
arma::vec params2vec_covar(const arma::mat& beta, int num_states, int Z_n_cols){
  arma::vec params_vec(num_states*(Z_n_cols)); params_vec.fill(0);
  int index = 0;
  for(int i = 0; i < beta.n_rows; i++){
    for(int j = 0; j < beta.n_cols; j++){
      params_vec(index) = beta(i,j);
      index++;
    }
  }
  return(params_vec);
}

// [[Rcpp::export]]
void get_emission_distribution(const arma::cube& Xt, const arma::cube& u_all,const arma::vec& m,
                               arma::cube& SampCov, arma::mat& mu_hat_next, arma::cube& Sigma_hat_next, arma::cube& Omega_hat_next){
  int N = Xt.n_slices;
  int P = Xt.n_cols;
  int TIME = Xt.n_rows;
  int num_states = m.n_elem;
  // int num_Bstates = accu(m);

  // Reset to 0
  SampCov.fill(0);
  mu_hat_next.fill(0);
  Sigma_hat_next.fill(0);
  Omega_hat_next.fill(0);

  for(int n = 0; n < N; n++){
    arma::mat Xtn = Xt.slice(n);
    arma::mat un = u_all.slice(n);

    for(int j = 0; j < num_states; j++){
      arma::vec uj(TIME); uj.fill(0);
      int id, K;
      if(j == 0){
        id = 0;
        K = m(0);
      }
      else{
        id = accu(m.subvec(0,j-1));
        K = accu(m.subvec(0,j));
      }
      for(int k = id; k < K; k++){
        uj += un.col(k);
      }
      //for(int k = id; k < K; k++){
        for(int p = 0; p < P; p++){
          mu_hat_next(j,p) += accu(uj%Xtn.col(p))/accu(uj)/N; // % is element-wise multiplication
          // mu_hat_next(k,p) += accu(un.col(k)%Xtn.col(p))/accu(un.col(k))/N; // % is element-wise multiplication
        }
      //}
    }
  }


  for(int n = 0; n < N; n++){
    arma::mat Xtn = Xt.slice(n);
    arma::mat un = u_all.slice(n);

    for(int j = 0; j < num_states; j++){
      arma::vec uj(TIME); uj.fill(0);
      int id, K;
      if(j == 0){
        id = 0;
        K = m(0);
      }
      else{
        id = accu(m.subvec(0,j-1));
        K = accu(m.subvec(0,j));
      }
      for(int k = id; k < K; k++){
        uj += un.col(k);
      }
      //for(int k = id; k < K; k++){
        SampCov.slice(j) += covmat_c(Xtn,mu_hat_next.row(j).t(),uj)/N;
        // SampCov.slice(k) += covmat_c(Xtn,mu_hat_next.row(k).t(),un.col(k))/N;
      //}
    }
  }
  Sigma_hat_next = SampCov;

}

// [[Rcpp::export]]
List forward_backward_hsmm_cpp(const arma::mat& Xt,
                               const arma::vec& m,
                         const arma::sp_mat& B,
                         const arma::mat& mu,
                         const arma::cube& root,
                         const arma::vec& delta,
                         bool chol=false){
  // B is transition probability matrix
  // mu is KxP matrix containing means, each row being the mean vector for that state
  // if chol = TRUE, root is PxPxK array, each slice being that state's square root of the covariance matrix
  // if chol = FALSE, root is PxPxK array, each slice being that state's covariance matrix
  int num_Bstates = accu(m);
  int num_states = m.n_elem;
  int TIME = Xt.n_rows;

  double fb_verysmallvalue = pow(4.940656,-142);

  // Get stationary distribution
  arma::mat ones(num_Bstates,num_Bstates); ones.fill(1);

  // Compute pdf at each time point
  arma::mat allprobs(TIME,num_Bstates); allprobs.fill(fb_verysmallvalue);
  for(int j = 0; j < num_states; j++){
    int id, K;
    if(j == 0){
      id = 0;
      K = m(0);
    }
    else{
      id = accu(m.subvec(0,j-1));
      K = accu(m.subvec(0,j));
    }
    allprobs.col(id) = dmvnrm_arma(Xt,mu.row(j),root.slice(j),FALSE,chol);

    for(int k = id+1; k < K; k++){
      allprobs.col(k) = allprobs.col(id);
    }
  }

  // Check if valid probability
  for(int time = 0; time < TIME; time++){
    if(all(allprobs.row(time) <= 1e-70)){
      allprobs.row(time).fill(1);
    }
    for(int j = 0; j < num_Bstates; j++){
      if(allprobs(time,j) < 0){
        allprobs(time,j) = -999;
        cout << "probability < 0" << endl;
      }
      if(!std::isfinite(allprobs(time,j))){
        allprobs(time,j) = -999;
        cout << "probability NaN" << endl;
      }
    }
  }

  // Forward variables
  arma::mat lalpha(num_Bstates,TIME); lalpha.fill(fb_verysmallvalue);
  arma::vec foo = delta%allprobs.row(0).t(); // element-wise multiplication
  double sumfoo = accu(foo) + fb_verysmallvalue;
  foo = foo/sumfoo;
  double lscale = std::log(sumfoo);
  lalpha.col(0) = logvec_c(foo) + lscale;
  for(int time = 1; time < TIME; time++){
    arma::sp_mat tempB(B.n_rows,B.n_cols);
    for(int i = 0; i < B.n_rows; i++){
      tempB.row(i) = (B.row(i)%allprobs.row(time));
    }
    foo = (foo.t()*tempB).t();
    sumfoo = accu(foo) + fb_verysmallvalue;
    lscale = lscale + std::log(sumfoo);
    foo = foo/accu(foo);
    lalpha.col(time) = logvec_c(foo) + lscale;
  }
  double llk_alpha = lscale;

  // Backward variables
  arma::mat lbeta(num_Bstates,TIME); lbeta.fill(fb_verysmallvalue);
  arma::vec foo2(num_Bstates); foo2.fill((double)(1/(double)num_Bstates));
  lscale = std::log(num_Bstates);
  for(int time = TIME-2; time >= 0; time--){
    foo2 = B*(allprobs.row(time+1)%foo2.t()).t();
    lbeta.col(time) = logvec_c(foo2) + lscale;
    double sumfoo2 = accu(foo2) + fb_verysmallvalue;
    foo2 = foo2/sumfoo2;
    lscale = lscale + std::log(sumfoo2);
  }
  double llk_beta = lscale;
  return(List::create(_["lalpha"] = lalpha,
                      _["lbeta"] = lbeta,
                      _["llk_alpha"] = llk_alpha,
                      _["llk_beta"] = llk_beta));
}

// [[Rcpp::export]]
long double pk_cpp(int r, double lambda, int tau = 1){
  // Shifted poisson pmf
  // lambda - rate parameter
  // tau - shift parameter
  return(R::dpois(r-tau,lambda,0));
  // return(pow(lambda,r-tau)*exp(-lambda)/factorial_cpp(r-tau));
}

// [[Rcpp::export]]
long double Fk_cpp(double r, double lambda, int tau = 1){
  // Shifted poisson cdf
  // lambda - rate parameter
  // tau - shift parameter
  return(R::ppois(r-tau,lambda,0,0)); // right tail for the cdf calculation
}

// [[Rcpp::export]]
long double ck_cpp(double r, double lambda, int tau = 1,int verbose = 0){
  if(r < tau){
    return(0);
  }
  if(std::isfinite(pk_cpp(r,lambda,tau)/(Fk_cpp(r-1,lambda,tau))) == 1){
    return(pk_cpp(r,lambda,tau)/(Fk_cpp(r-1,lambda,tau)));
  }
  if(Fk_cpp(r,lambda,tau) == 0){
    return(1);
  }
  if(verbose == 1){
    cout << "Error in ck_cpp!" << endl;
  }
  return(-1);
}


// [[Rcpp::export]]
long double grad_ck_cpp(double r, double lambda, int tau = 1,int verbose = 0){
  if(r < tau){
    return(0);
  }
  if(std::isfinite(pk_cpp(r,lambda,tau)/(Fk_cpp(r-1,lambda,tau))) == 1){
    return((-pk_cpp(r,lambda,tau)*pk_cpp(r-1,lambda,tau)+Fk_cpp(r-1,lambda,tau)*(pk_cpp(r-1,lambda,tau)-pk_cpp(r,lambda,tau)))/pow(Fk_cpp(r-1,lambda,tau),2));
  }
  if(Fk_cpp(r,lambda,tau) == 0){
    return(0);
  }
  if(verbose == 1){
    cout << "Error in ck_cpp!" << endl;
  }
  return(-1);
}

// [[Rcpp::export]]
arma::mat create_Bij_cpp(double lambda, double aij, int mi, int mj, int tau = 1){
  arma::mat Bij(mi,mj); Bij.fill(0);
  arma::vec Bijvec(mi); Bijvec.fill(0);
  for(int j = 0; j < mi; j++){
    Bijvec(j) = aij*ck_cpp(j+1,lambda,tau);
  }
  Bij.col(0) = Bijvec;
  return(Bij);
}

// [[Rcpp::export]]
arma::mat create_gradBij_cpp(double lambda, int mi, int mj, int tau = 1){
  arma::mat gradBij(mi,mj); gradBij.fill(0);
  arma::vec gradBijvec(mi); gradBijvec.fill(0);
  for(int j = 0; j < mi; j++){
    gradBijvec(j) = pow(ck_cpp(j+1,lambda,tau),-1)*grad_ck_cpp(j+1,lambda,tau);
  }
  gradBij.col(0) = gradBijvec;
  return(gradBij);
}

// [[Rcpp::export]]
arma::mat create_Bii_cpp(double lambda, double mi, int tau=1,int verbose = 0){
  arma::vec diagB(mi-1); diagB.fill(0);
  for(int i = 0; i < (mi-1); i++){
    diagB(i) = 1-ck_cpp(i+1,lambda,tau, verbose);
  }
  arma::mat Bii(mi,mi); Bii.fill(0);
  Bii.submat(0,0,mi-2,0).fill(0);
  Bii.submat(0,1,mi-2,mi-1) = diagmat(diagB);

  arma::vec Bii_lastrow(mi); Bii_lastrow.fill(0);
  Bii_lastrow(mi-1) = 1-ck_cpp(mi,lambda,tau);
  Bii.row(mi-1) = Bii_lastrow.t();
  return(Bii);
}

// [[Rcpp::export]]
arma::mat create_gradBii_cpp(double lambda, double mi, int tau=1,int verbose = 0){
  arma::vec diaggradB(mi-1); diaggradB.fill(0);
  for(int i = 0; i < (mi-1); i++){
    diaggradB(i) = -pow(1-ck_cpp(i+1,lambda,tau, verbose),-1)*grad_ck_cpp(i+1,lambda,tau);
  }
  arma::mat gradBii(mi,mi); gradBii.fill(0);
  gradBii.submat(0,0,mi-2,0).fill(0);
  gradBii.submat(0,1,mi-2,mi-1) = diagmat(diaggradB);

  arma::vec gradBii_lastrow(mi); gradBii_lastrow.fill(0);
  gradBii_lastrow(mi-1) = -pow(1-ck_cpp(mi,lambda,tau),-1)*grad_ck_cpp(mi,lambda,tau);
  gradBii.row(mi-1) = gradBii_lastrow.t();
  return(gradBii);
}

// [[Rcpp::export]]
arma::sp_mat create_B_cpp(arma::vec lambda, arma::vec m, arma::mat A,int tau = 1,int verbose = 0){
  int num_states = m.n_elem;
  arma::mat B(accu(m),accu(m)); B.fill(0);

  // First "state aggregate"
  B.submat(0,0,m(0)-1,m(0)-1) = create_Bii_cpp(lambda(0),m(0), tau,verbose);
  for(int i = 1; i < num_states; i++){
    B.submat(0,accu(m.subvec(0,i-1)),
             m(0)-1,accu(m.subvec(0,i))-1) = create_Bij_cpp(lambda(0),A(0,i),m(0),m(i));
  }

  // The other state aggregate
  for(int i = 1; i < num_states; i++){
    for(int j = 0; j < num_states; j++){
      // Special treatment for the first state aggregate
      if(j==0){
        B.submat(accu(m.subvec(0,i-1)),0,
                 accu(m.subvec(0,i))-1,m(j)-1) = create_Bij_cpp(lambda(i),A(i,j),m(i),m(j));
      }
      // Diagonal
      if(i == j){
        B.submat(accu(m.subvec(0,i-1)),accu(m.subvec(0,i-1)),
                 accu(m.subvec(0,i))-1,accu(m.subvec(0,i))-1) = create_Bii_cpp(lambda(i),m(i));
      }
      // Left and right of the diagonal
      if(i != j && j != 0){
        B.submat(accu(m.subvec(0,i-1)),accu(m.subvec(0,j-1)),
                 accu(m.subvec(0,i))-1,accu(m.subvec(0,j))-1) = create_Bij_cpp(lambda(i),A(i,j),m(i),m(j));
      }
    }
  }
  return(arma::sp_mat(B));
}

// [[Rcpp::export]]
arma::mat create_gradB_cpp(arma::vec lambda, arma::vec m, int tau = 1,int verbose = 0){
  int num_states = m.n_elem;
  arma::mat gradB(accu(m),accu(m)); gradB.fill(0);

  // First "state aggregate"
  gradB.submat(0,0,m(0)-1,m(0)-1) = create_gradBii_cpp(lambda(0),m(0), tau,verbose);
  for(int i = 1; i < num_states; i++){
    gradB.submat(0,accu(m.subvec(0,i-1)),
             m(0)-1,accu(m.subvec(0,i))-1) = create_gradBij_cpp(lambda(0),m(0),m(i));
  }

  // The other state aggregate
  for(int i = 1; i < num_states; i++){
    for(int j = 0; j < num_states; j++){
      // Special treatment for the first state aggregate
      if(j==0){
        gradB.submat(accu(m.subvec(0,i-1)),0,
                 accu(m.subvec(0,i))-1,m(j)-1) = create_gradBij_cpp(lambda(i),m(i),m(j));
      }
      // Diagonal
      if(i == j){
        gradB.submat(accu(m.subvec(0,i-1)),accu(m.subvec(0,i-1)),
                 accu(m.subvec(0,i))-1,accu(m.subvec(0,i))-1) = create_gradBii_cpp(lambda(i),m(i));
      }
      // Left and right of the diagonal
      if(i != j && j != 0){
        gradB.submat(accu(m.subvec(0,i-1)),accu(m.subvec(0,j-1)),
                 accu(m.subvec(0,i))-1,accu(m.subvec(0,j))-1) = create_gradBij_cpp(lambda(i),m(i),m(j));
      }
    }
  }
  return(gradB);
}

// [[Rcpp::export]]
arma::vec get_lambdan(const arma::vec& Zn, const arma::mat& beta){
  int num_states = beta.n_rows;
  arma::vec lambda(num_states); lambda.fill(0);
  for(int i = 0; i < num_states; i++){
    double loglambda = 0;
    for(int j = 0; j < Zn.n_elem; j++){
      loglambda += beta(i,j)*Zn(j);
    }
    lambda(i) = exp(loglambda);
  }
  return lambda;
}

// [[Rcpp::export]]
double term2_covar_cpp(arma::vec params_vec,
                       const arma::cube& Amat,
                       arma::mat Z, arma::vec m, const arma::cube& sum_t_vjk,
                       int verbose = 0){
  int N = sum_t_vjk.n_slices;
  int num_states = m.n_elem;
  // int num_Bstates = accu(m);
  int tau = 1; // Shifted away from 0 since we are forcing a non-zero duration

  // mat ones(num_Bstates,num_Bstates); ones.fill(1);
  arma::mat betahat(num_states,Z.n_cols); betahat.fill(0);
  vec2params_covar(params_vec,num_states,Z.n_cols,betahat);

  double term2 = 0;
  for(int n = 0; n < N; n++){
    arma::vec lambda = get_lambdan(Z.row(n).t(), betahat);

    // Get tpm
    arma::sp_mat B = create_B_cpp(lambda,m,Amat.slice(n),tau,verbose);

    // Term 2
    for(int j = 0; j < B.n_rows; j++){
      for(int k = 0; k < B.n_cols; k++){
        //if(B_id(i,j) == 1 && sum_t_vjk(i,j,n) > 0){
        if(B(j,k) > 0){
          term2 += std::log(B(j,k))*sum_t_vjk(j,k,n);
        }
      }
    }
  }
  return(-(term2));
}

// [[Rcpp::export]]
arma::mat grad_term2_cpp(arma::vec params_vec,
                   const arma::cube& Amat,
                      const arma::mat& Z,
                      const arma::vec& m,
                      const arma::cube sum_t_vjk){
  // Gradient of term2_covar_cpp
  //int tau = 1;
  int N = sum_t_vjk.n_slices;
  int num_states = m.n_elem;
  arma::mat betahat(num_states,Z.n_cols); betahat.fill(0);
  vec2params_covar(params_vec,num_states,Z.n_cols,betahat);
  arma::mat gradbeta(num_states,Z.n_cols); gradbeta.fill(0);

  for(int n = 0; n < N; n++){
    arma::vec lambda = get_lambdan(Z.row(n).t(), betahat);
    arma::mat gradB = create_gradB_cpp(lambda,m);
    arma::mat sum_t_vjkn = sum_t_vjk.slice(n);
    for(int j = 0; j < num_states; j++){
      for(int k = 0; k < Z.n_cols; k++){
        if(j == 0){
          gradbeta(j,k) += accu(sum_t_vjkn.submat(0,0,m[1]-1,sum_t_vjkn.n_cols-1)%gradB.submat(0,0,m[1]-1,sum_t_vjkn.n_cols-1)*Z(n,k)*lambda(j));
        }
        if(j > 0){
          gradbeta(j,k) += accu(sum_t_vjkn.submat(accu(m.subvec(0,j-1)),0,accu(m.subvec(0,j))-1,sum_t_vjkn.n_cols-1)%gradB.submat(accu(m.subvec(0,j-1)),0,accu(m.subvec(0,j))-1,sum_t_vjkn.n_cols-1)*Z(n,k)*lambda(j));
        }
      }
    }
  }
  return(-gradbeta); // nevative because we're maximizing, not minimizing, and optim minimizes
}

// [[Rcpp::export]]
arma::vec optim_covar_rcpp(const arma::vec& params_vec, const arma::cube& Amat,
                           const arma::mat& Z, const arma::vec& m, const arma::cube& sum_t_vjk){
  // Extract R's optim function

  Rcpp::Environment stats("package:stats");
  Rcpp::Function optim = stats["optim"];

  Rcpp::List results = optim(Rcpp::_["par"]    = params_vec,
                             Rcpp::_["fn"]     = Rcpp::InternalFunction(&term2_covar_cpp),
                             Rcpp::_["method"] = "Nelder-Mead",
                             Rcpp::_["hessian"] = 0,
                             Rcpp::_["Amat"] = Amat,
                             //Rcpp::_["tau"] = tau,
                             Rcpp::_["Z"] = Z,
                             Rcpp::_["m"] = m,
                             Rcpp::_["sum_t_vjk"] = sum_t_vjk,
                             Rcpp::_["verbose"] = 0);

  /*
  Rcpp::List results = optim(Rcpp::_["par"]    = params_vec,
                             Rcpp::_["fn"]     = Rcpp::InternalFunction(&term2_covar_cpp),
                             Rcpp::_["gr"]     = Rcpp::InternalFunction(&grad_term2_cpp),
                             Rcpp::_["method"] = "BFGS",
                             Rcpp::_["hessian"] = 0,
                             Rcpp::_["Amat"] = Amat,
                             //Rcpp::_["tau"] = tau,
                             Rcpp::_["Z"] = Z,
                             Rcpp::_["m"] = m,
                             Rcpp::_["sum_t_vjk"] = sum_t_vjk,
                             Rcpp::_["verbose"] = 0);
   */

  arma::vec out = Rcpp::as<arma::vec>(results[0]);

  return out;
}

// [[Rcpp::export]]
void em_estep_covar(const arma::cube& Xt, const arma::mat& Z, const arma::vec& m,
                    const arma::mat& beta, const arma::cube& A, const arma::mat& mu, const arma::cube& Sigma,const arma::mat& delta, int iter,
              arma::cube& lalpha, arma::cube& lbeta, arma::mat& llk, arma::cube& sum_t_vjk, arma::cube& u_all){
  // E step
  //////
  int TIME = Xt.n_rows;
  int N = Xt.n_slices;
  int num_Bstates = accu(m);
  int num_states = m.size();

  // Do the choleskys in advance
  arma::cube root_all(Sigma.n_rows,Sigma.n_cols,num_states);
  root_all.slice(0) = arma::trans(arma::inv(trimatu(arma::chol(Sigma.slice(0)))));
  for(int i = 1; i < num_states; i++){
    root_all.slice(i) = arma::trans(arma::inv(trimatu(arma::chol(Sigma.slice(i)))));
  }

  // Now do fb algorithm for each subject
  for(int n = 0; n < N; n++){
    arma::mat Xtn = Xt.slice(n);
    arma::sp_mat B = create_B_cpp(get_lambdan(Z.row(n).t(), beta),
                            m,
                            A.slice(n));

    // Forward-backward algorithm
    List fb = forward_backward_hsmm_cpp(Xtn,m,B,mu,root_all,delta.row(n).t(),TRUE);

    // Calculate likelihood
    arma::mat lalpha_n = fb[0];
    arma::mat lbeta_n = fb[1];

    lalpha.slice(n) = lalpha_n;
    lbeta.slice(n) = lbeta_n;
    arma::mat lprobs(TIME,num_Bstates); lprobs.fill(0);
    for(int i = 0; i < num_states; i++){
      int index = 0;
      if(i > 0){
        index = accu(m.subvec(0,i-1));
      }
      else{
        index = 0;
      }
      lprobs.col(index) = dmvnrm_arma(Xtn,
                 mu.row(i),
                 root_all.slice(i),
                 TRUE,
                 TRUE);
      for(int j = 1; j < m(i); j++){
        lprobs.col(index + j) = lprobs.col(index);
      }
    }
    //cout << "lprobs:" << endl;
    //lprobs.print();

    // Get transition probabilities and the "u" vector (estimate of state sequence)
    arma::mat u_mat(TIME,num_Bstates); u_mat.fill(0);
    llk(n,iter) = fb[2];

    for(int j = 0; j < num_Bstates; j++){
      for(int k = 0; k < num_Bstates; k++){
        //if(B(j,k) > pow(1,-16)){
          arma::vec lalpha_temp = lalpha_n.row(j).subvec(0,TIME-2).t();
          arma::vec lbeta_temp = lbeta_n.row(k).subvec(1,TIME-1).t();
          arma::vec lprobs_temp = lprobs.col(k).subvec(1,TIME-1);
          arma::vec logBjk_temp(TIME-1); logBjk_temp.fill(std::log(B(j,k)));
          sum_t_vjk(j,k,n) = accu(exp(lalpha_temp+lbeta_temp+logBjk_temp+lprobs_temp-llk(n,iter)));
      }
      // Now get the u vector
      rowvec u = exp(lalpha_n.row(j) + lbeta_n.row(j) - llk(n,iter));
      u_mat.col(j) = u.t();
    }
    u_all.slice(n) = u_mat;
  }
}

// [[Rcpp::export]]
arma::cube get_Ahat(const arma::cube& sum_t_vjk,
               const arma::vec& m,
               int num_states){
  arma::cube Ahat(num_states,num_states,sum_t_vjk.n_slices); Ahat.fill(0);
  for(int n = 0; n < sum_t_vjk.n_slices; n++){
    arma::mat vbar(num_states,num_states); vbar.fill(0); // Stores sums of sum_t_vjk
    arma::mat sum_t_vjkn = sum_t_vjk.slice(n);
    for(int s = 0; s < num_states; s++){
      // Upper triangle
      if(s != num_states-1){ // No need to do calculations if working on last row
        for(int r = (s+1); r < num_states; r++){
          if(s == 0){
            vbar(s,r) += accu(sum_t_vjkn.submat(0,
                              accu(m.subvec(0,r-1)),
                              m(0)-1,
                              accu(m.subvec(0,r-1))));
          }
          else{
            vbar(s,r) += accu(sum_t_vjkn.submat(accu(m.subvec(0,s-1)),
                              accu(m.subvec(0,r-1)),
                              accu(m.subvec(0,s))-1,
                              accu(m.subvec(0,r-1))));
          }
        }
      }

      // Lower triangle
      if(s != 0){// No need to do calculations if working on first row
        for(int r = 0; r < s; r++){
          if(r == 0){
            vbar(s,r) += accu(sum_t_vjkn.submat(accu(m.subvec(0,s-1)),
                              0,
                              accu(m.subvec(0,s))-1,
                              0));
          }
          else{
            vbar(s,r) += accu(sum_t_vjkn.submat(accu(m.subvec(0,s-1)),
                              accu(m.subvec(0,r-1)),
                              accu(m.subvec(0,s))-1,
                              accu(m.subvec(0,r-1))));
          }
        }
      }
    }

    for(int s = 0; s < num_states; s++){
      Ahat.slice(n).row(s) = vbar.row(s)/accu(vbar.row(s));
    }
  }

  return(Ahat);
}

// [[Rcpp::export]]
void em_mstep_covar(const arma::cube& Xt, const arma::mat& Z, const arma::vec& m,
                    const arma::cube& sum_t_vjk, const arma::cube& u_all,
                    const arma::mat& beta_hat, arma::mat& beta_hat_next, arma::cube& A_hat_next, arma::mat& delta_hat_next,
                    arma::mat& mu_hat_next, arma::cube& SampCov, arma::cube& Sigma_hat_next, arma::cube& Omega_hat_next, bool refit){
  // Maximize conditional likelihood
  int num_states = m.n_elem;

  if(!refit){
    // Get TPM
    A_hat_next = get_Ahat(sum_t_vjk,m,num_states);

    // Get stationary distribution
    for(int n = 0; n < u_all.n_slices; n++){
      arma::mat un = u_all.slice(n);
      delta_hat_next.row(n) = un.row(0)/accu(un.row(0));
      for(int j = 0; j < delta_hat_next.n_cols; j++){
        if(delta_hat_next(n,j) <= 0){
          delta_hat_next(n,j) = pow(4.940656,-142);
        }
      }
    }
  }


  // Get dwell-time distribution parameters
  arma::vec params_vec = params2vec_covar(beta_hat,num_states,Z.n_cols);
  arma::vec params_hat = optim_covar_rcpp(params_vec,
                                    A_hat_next,Z,m,
                                    sum_t_vjk);

  vec2params_covar(params_hat,num_states,Z.n_cols,beta_hat_next);

  // Get emission distribution parameters
  // Mean vector, covariance matrix, precision matrix
  if(!refit){
    get_emission_distribution(Xt, u_all,m,
                              SampCov, mu_hat_next, Sigma_hat_next, Omega_hat_next);
  }
}

// [[Rcpp::export]]
List hsmm_cpp(const arma::cube& Xt,
                            const arma::mat& Z,
                            const arma::cube& A_init,
                            const arma::mat& mu_init,
                            const arma::cube& Sigma_init,
                            const arma::mat& beta_init,
                            const arma::vec& m,
                            int maxiter = 10000, double tol = 1e-4, bool verbose = 1,
                            bool refit = 0){
  int TIME = Xt.n_rows;
  int N = Xt.n_slices;
  int P = Xt.n_cols;
  int num_states = A_init.n_rows;
  int num_Bstates = accu(m);

  arma::vec lambda(beta_init.n_rows); lambda.fill(0);
  arma::cube sum_t_vjk(num_Bstates,num_Bstates,N); sum_t_vjk.fill(0);

  // Set ''current'' parameters to initialized parameters
  arma::mat mu_hat = mu_init;
  arma::cube Sigma_hat = Sigma_init;
  arma::mat beta_hat = beta_init;
  // Transform A so that diagonals are 0 and off-diagonals add up to 1, row-wise
  arma::cube A_hat = A_init;
  for(int n = 0; n < A_hat.n_slices; n++){
    A_hat.slice(n).diag().fill(0);
    A_hat.slice(n) = standardize_rows_cpp(A_hat.slice(n));
  }

  // Initialize parameters
  arma::mat mu_hat_next(num_states,P); mu_hat_next.fill(0);
  arma::cube Sigma_hat_next(P,P,num_states); Sigma_hat_next.fill(0);
  arma::cube Omega_hat_next(P,P,num_states); Omega_hat_next.fill(0);
  arma::mat beta_hat_next(num_states,Z.n_cols); beta_hat_next.fill(0);
  arma::cube A_hat_next(num_states,num_states,N); A_hat_next.fill(0);
  arma::mat delta_hat(N,num_Bstates); delta_hat.fill(1/(double)num_Bstates);
  arma::mat delta_hat_next(N,num_Bstates); delta_hat_next.fill(0);
  arma::cube SampCov(P,P,num_states); SampCov.fill(0);
  arma::mat llk(N,maxiter); llk.fill(-999);
  arma::cube u_all(TIME,num_Bstates,N); u_all.fill(0);
  arma::cube lalpha(num_Bstates,TIME,N); lalpha.fill(0);
  arma::cube lbeta(num_Bstates,TIME,N); lbeta.fill(0);
  arma::mat Omega_hat(P,P); Omega_hat.fill(0);
  for(int iter = 0; iter < maxiter; iter++){
    if(verbose == 1){
      cout << "Iteration " << iter+1 << " of " << maxiter << "..." << endl;
    }

    // Reset parameters
    if(refit){
      beta_hat_next.fill(0);
      // Hold these parameters constant
      A_hat_next = A_hat;
      delta_hat_next = delta_hat;
      mu_hat_next = mu_hat;
      Sigma_hat_next = Sigma_hat;
    }
    else{
      A_hat_next.fill(0);
      mu_hat_next.fill(0);
      Sigma_hat_next.fill(0);
      SampCov.fill(0);
      beta_hat_next.fill(0);
    }

    // E step
    /////
    // cout << "E step" << endl;
    em_estep_covar(Xt, Z,m,
                   beta_hat, A_hat, mu_hat, Sigma_hat, delta_hat,
                   iter,
                   lalpha, lbeta, llk, sum_t_vjk, u_all);

    // M step
    /////
    // cout << "M step" << endl;
    em_mstep_covar(Xt, Z, m,
                   sum_t_vjk, u_all,
                   beta_hat, beta_hat_next, A_hat_next, delta_hat_next,
                   mu_hat_next, SampCov, Sigma_hat_next, Omega_hat_next,
                   refit);

    if(refit){
      SampCov = Sigma_init;
      Sigma_hat_next = Sigma_init;
      mu_hat_next = mu_init;
    }

    A_hat = A_hat_next;
    delta_hat = delta_hat_next;
    beta_hat = beta_hat_next;
    mu_hat = mu_hat_next;
    Sigma_hat = Sigma_hat_next;

    // Determine convergence
    if(iter > 0){
      if(verbose == 1){
        cout << "Likelihood this iteration = " << accu(llk.col(iter)) << endl;
        cout << "Likelihood previous iteration = " << accu(llk.col(iter-1)) << endl;
        cout << "Relative Change in likelihood = " << (accu(llk.col(iter)-llk.col(iter-1)))/fabs(accu(llk.col(iter-1))) << endl;
      }
      bool converged = fabs((accu(llk.col(iter)-llk.col(iter-1)))/accu(llk.col(iter-1))) < tol;
      if(converged){
        cout << "Converged after " << iter+1 << " iterations" << endl;
          return(List::create(_["mu.hat"] = mu_hat_next,
                              _["Sigma.hat"] = Sigma_hat_next,
                              _["SampCov"] = SampCov,
                              _["A.hat"] = A_hat_next,
                              _["beta.hat"] = beta_hat_next,
                              _["delta.hat"] = delta_hat_next,
                              _["llk"] = llk.submat(0,0,N-1,iter),
                              _["u.all"] = u_all,
                              _["sum_t_vjk"] = sum_t_vjk));
      }
    }

  }
  cout << "No convergence after " << maxiter << " iterations" << endl;
  return(List::create(_["mu.hat"] = mu_hat_next,
                      _["Sigma.hat"] = Sigma_hat_next,
                      _["SampCov"] = SampCov,
                      _["A.hat"] = A_hat_next,
                      _["beta.hat"] = beta_hat_next,
                      _["llk"] = llk,
                      _["u.all"] = u_all));
}
