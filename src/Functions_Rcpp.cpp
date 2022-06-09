#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Inverse of a symmetric matrix
// [[Rcpp::export]]
arma::mat inversa(arma::mat M) {
  return (arma::inv(M));
} 


// Quatities to estimate beta
// [[Rcpp::export]]
List ciclobeta(arma::mat X, arma::vec phi, arma::vec SAEMu, arma::vec SAEMuyi, arma::mat SAEMuzi){
  int n1 = SAEMu.n_elem; int q = X.n_cols; int p = phi.n_elem; 
  arma::vec A0(q,fill::zeros); arma::mat Ai(q,q,fill::zeros); arma::vec Bi(q,fill::zeros);
  arma::uvec index(p,fill::zeros);
  
  for (int i=0; i<n1; i++){
    std::iota(index.begin(),index.end(),i);
    std::reverse(index.begin(),index.end());
    A0 = (X.row(p+i)).t() - (X.rows(index)).t()*phi;
    Ai += SAEMu(i)*(A0*A0.t());
    Bi += (SAEMuyi(i) - as_scalar(phi.t()*SAEMuzi.col(i)))*A0;
  }
  
  List ciclos;
  ciclos["Ai"] = Ai; ciclos["Bi"] = Bi;
  return ciclos;
}


// Compute the mean of Y|Y_p, u, theta
// [[Rcpp::export]]
List ComputeMean(arma::vec yobs, arma::mat Xmatrix, arma::vec beta, arma::vec phi){
  int nj = yobs.n_elem; int pj = phi.n_elem; 
  arma::mat diagm((pj-1),(pj-1),fill::eye); arma::vec ov(pj-1,fill::zeros); 
  arma::mat Phimatrix = join_vert(phi.t(),join_horiz(diagm,ov));
  arma::mat Phi(nj-pj,pj,fill::zeros);
  
  arma::uvec index(pj,fill::zeros);
  std::iota(index.begin(),index.end(),0);
  std::reverse(index.begin(),index.end());
  arma::vec meanp =  yobs.elem(index) - Xmatrix.rows(index)*beta;
  
  arma::mat pom;
  for (int i=0; i<(nj-pj); i++){
    pom = powmat(Phimatrix,(i+1));
    Phi.row(i) = (pom).row(0);
  }
  
  arma::uvec index2(nj-pj,fill::zeros);
  std::iota(index2.begin(),index2.end(),pj);
  arma::vec meanvc = Xmatrix.rows(index2)*beta + Phi*meanp;
  List mediacond;
  mediacond["mvc"] = meanvc;
  mediacond["Pphi"] = Phi;
  return(mediacond);
}


// Compute the variance of Y|Y_p, u, theta
// [[Rcpp::export]]
arma::mat ComputeVar(double sigma2, arma::mat Pphi, arma::vec phi, arma::vec SAEMu){
  int mj = SAEMu.n_elem; int pj = phi.n_elem;
  arma::mat varvc(mj,mj,fill::zeros); arma::rowvec ab(pj, fill::zeros); ab(0) = 1; 
  arma::mat Phimatrix = join_vert(ab, Pphi); double sum1;
  
  for (int i=0; i<mj; i++){
    for (int j=i; j<mj; j++){
      sum1 = 0;
      for (int k=0; k<(i+1); k++){
        sum1 += Phimatrix((i-k),0)*Phimatrix((j-k),0)/SAEMu(k);
      }
      varvc(i,j) = varvc(j,i) = sum1*sigma2;
    }
  }
  return varvc;
}


// Score vector of the complete log-likelihood function
// [[Rcpp::export]]
arma::vec Gradient(arma::mat x, arma::vec saemU, arma::vec saemY, arma::vec phi, arma::vec beta, 
                   double sigma2, double nu, bool fixnu){
  int n = saemY.n_elem; int p = phi.n_elem; int n1 = n-p; int q = beta.n_elem;
  double sum1 = 0; arma::vec sum21(p,fill::zeros); arma::mat sum22(p,p,fill::zeros);  
  arma::vec sum31(q,fill::zeros); arma::mat sum32(q,q,fill::zeros);
  arma::vec ytilde = saemY - x*beta;
  arma::uvec index(p,fill::zeros); arma::vec score;

  for (int i=0; i<n1; i++){
    std::iota(index.begin(),index.end(),i);
    std::reverse(index.begin(),index.end());
    sum1 += saemU(i)*pow(ytilde(p+i) - as_scalar((ytilde.elem(index)).t()*phi) ,2);
    sum21 += saemU(i)*ytilde(p+i)*ytilde.elem(index);
    sum22 += saemU(i)*(ytilde.elem(index)*(ytilde.elem(index)).t());
    sum31 += saemU(i)*(saemY(p+i) - as_scalar((saemY.elem(index)).t()*phi))*((x.row(p+i)).t() - (x.rows(index)).t()*phi);
    sum32 += saemU(i)*((x.row(p+i)).t() - (x.rows(index)).t()*phi)*((x.row(p+i)).t() - (x.rows(index)).t()*phi).t();
  }
  // First derivatives
  if (fixnu==false){
    arma::vec dbeta; arma::vec dphi; arma::vec dsigma; arma::vec dnu;
    
    dbeta = (1/sigma2)*(sum31 - sum32*beta);
    dphi = (1/sigma2)*(sum21 - sum22*phi);
    dsigma = (-0.5*n1/sigma2) + (0.5/pow(sigma2,2))*sum1;
    dnu = 0.5*n1*(log(0.5*nu) + 1 - R::digamma(0.5*nu)) + 0.5*(sum(log(saemU)) - sum(saemU));
    
    score = join_vert(dbeta, dsigma, dphi, dnu);
  } else {
    arma::vec dbeta; arma::vec dphi; arma::vec dsigma;
    
    dbeta = (1/sigma2)*(sum31 - sum32*beta);
    dphi = (1/sigma2)*(sum21 - sum22*phi);
    dsigma = (-0.5*n1/sigma2) + (0.5/pow(sigma2,2))*sum1;
    
    score = join_vert(dbeta, dsigma, dphi);
  }
  return score;
}


// Expectation of Score vector of the complete log-likelihood function
// [[Rcpp::export]]
arma::vec GradientExp(double saemLU, arma::vec saemU, double uy2, arma::vec uyw, arma::mat uw2, arma::mat A, 
                      arma::vec B, arma::vec phi, arma::vec beta, double sigma2, double nu, bool fixnu){
  int n1 = saemU.n_elem; arma::vec score;
  
  if (fixnu==false){
    arma::vec dbeta; arma::vec dphi; arma::vec dsigma; arma::vec dnu;
  
    dbeta = (1/sigma2)*(B - A*beta);
    dphi = (1/sigma2)*(uyw - uw2*phi);
    dsigma = (-0.5*n1/sigma2) + (0.5/pow(sigma2,2))*(uy2 - phi.t()*uyw - uyw.t()*phi + phi.t()*uw2*phi);
    dnu = 0.5*n1*(log(0.5*nu) + 1 - R::digamma(0.5*nu)) + 0.5*(saemLU - sum(saemU));

    score = join_vert(dbeta, dsigma, dphi, dnu);
  } else {
    arma::vec dbeta; arma::vec dphi; arma::vec dsigma;
    
    dbeta = (1/sigma2)*(B - A*beta);
    dphi = (1/sigma2)*(uyw - uw2*phi);
    dsigma = (-0.5*n1/sigma2) + (0.5/pow(sigma2,2))*(uy2 - phi.t()*uyw - uyw.t()*phi + phi.t()*uw2*phi);
    
    score = join_vert(dbeta, dsigma, dphi);
  }
  return score;
}


// Expectation of Hessian matrix of the complete log-likelihood function
// [[Rcpp::export]]
arma::mat HessianExp(arma::vec saemU, arma::mat saemUZi, arma::vec saemUYi, double uy2, arma::vec uyw, arma::mat uw2, 
                     arma::mat A, arma::vec B, arma::vec media, arma::mat X, arma::vec phi, arma::vec beta, 
                     double sigma2, double nu, bool fixnu){
  int n1 = saemU.n_elem; int p = phi.n_elem; int q = beta.n_elem;
  arma::mat sum22(q,p,fill::zeros); double b0 = 0;
  arma::uvec index(p,fill::zeros); arma::mat hessian;
  
  for (int i=0; i<n1; i++){
    std::iota(index.begin(),index.end(),i);
    std::reverse(index.begin(),index.end());
    b0 = as_scalar((saemUZi.col(i)).t()*phi) + saemU(i)*(media(p+i) - as_scalar((media.elem(index)).t()*phi)) - saemUYi(i);
    sum22 += b0*(X.rows(index)).t() - ((X.row(p+i)).t() - (X.rows(index)).t()*phi)*(saemUZi.col(i) - saemU(i)*media.elem(index)).t();
  }

  // Second derivatives
  if (fixnu == false){
    arma::mat d2beta; arma::mat dbetadphi; arma::vec dbetadsigma; arma::vec dbetadnu(q,fill::zeros); 
    arma::mat d2phi; arma::vec dphidsigma; arma::vec dphidnu(p,fill::zeros); arma::vec d2sigma;
    arma::vec dsigmadnu(1,fill::zeros); arma::vec d2nu;
    
    d2beta = (-1/sigma2)*A;
    dbetadphi = (1/sigma2)*sum22;
    dbetadsigma = (-1/pow(sigma2,2))*(B - A*beta);
    
    d2phi = (-1/sigma2)*uw2;
    dphidsigma = (-1/pow(sigma2,2))*(uyw - uw2*phi);
    
    d2sigma = (0.5*n1/pow(sigma2,2)) - (1/pow(sigma2,3))*(uy2 - phi.t()*uyw - uyw.t()*phi + phi.t()*uw2*phi);
    d2nu = (0.5*n1)*(1/nu - 0.5*R::trigamma(0.5*nu));
    
    hessian = join_vert(join_horiz(d2beta, dbetadsigma, dbetadphi, dbetadnu),
                                  join_horiz(dbetadsigma.t(), d2sigma, dphidsigma.t(), dsigmadnu),
                                  join_horiz(dbetadphi.t(), dphidsigma, d2phi, dphidnu),
                                  join_horiz(dbetadnu.t(), dsigmadnu, dphidnu.t(), d2nu));
  } else {
    arma::mat d2beta; arma::mat dbetadphi; arma::vec dbetadsigma;
    arma::mat d2phi; arma::vec dphidsigma; arma::vec d2sigma;
    
    d2beta = (-1/sigma2)*A;
    dbetadphi = (1/sigma2)*sum22;
    dbetadsigma = (-1/pow(sigma2,2))*(B - A*beta);
    d2phi = (-1/sigma2)*uw2;
    dphidsigma = (-1/pow(sigma2,2))*(uyw - uw2*phi);
    d2sigma = (0.5*n1/pow(sigma2,2)) - (1/pow(sigma2,3))*(uy2 - phi.t()*uyw - uyw.t()*phi + phi.t()*uw2*phi);

    hessian = join_vert(join_horiz(d2beta, dbetadsigma, dbetadphi),
                                  join_horiz(dbetadsigma.t(), d2sigma, dphidsigma.t()),
                                  join_horiz(dbetadphi.t(), dphidsigma, d2phi));
  }
  return hessian;
}
