// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// inversa
arma::mat inversa(arma::mat M);
RcppExport SEXP _ARCensReg_inversa(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(inversa(M));
    return rcpp_result_gen;
END_RCPP
}
// ciclobeta
List ciclobeta(arma::mat X, arma::vec phi, arma::vec SAEMu, arma::vec SAEMuyi, arma::mat SAEMuzi);
RcppExport SEXP _ARCensReg_ciclobeta(SEXP XSEXP, SEXP phiSEXP, SEXP SAEMuSEXP, SEXP SAEMuyiSEXP, SEXP SAEMuziSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type SAEMu(SAEMuSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type SAEMuyi(SAEMuyiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SAEMuzi(SAEMuziSEXP);
    rcpp_result_gen = Rcpp::wrap(ciclobeta(X, phi, SAEMu, SAEMuyi, SAEMuzi));
    return rcpp_result_gen;
END_RCPP
}
// ComputeMean
List ComputeMean(arma::vec yobs, arma::mat Xmatrix, arma::vec beta, arma::vec phi);
RcppExport SEXP _ARCensReg_ComputeMean(SEXP yobsSEXP, SEXP XmatrixSEXP, SEXP betaSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type yobs(yobsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xmatrix(XmatrixSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phi(phiSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeMean(yobs, Xmatrix, beta, phi));
    return rcpp_result_gen;
END_RCPP
}
// ComputeVar
arma::mat ComputeVar(double sigma2, arma::mat Pphi, arma::vec phi, arma::vec SAEMu);
RcppExport SEXP _ARCensReg_ComputeVar(SEXP sigma2SEXP, SEXP PphiSEXP, SEXP phiSEXP, SEXP SAEMuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Pphi(PphiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type SAEMu(SAEMuSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeVar(sigma2, Pphi, phi, SAEMu));
    return rcpp_result_gen;
END_RCPP
}
// Gradient
arma::vec Gradient(arma::mat x, arma::vec saemU, arma::vec saemY, arma::vec phi, arma::vec beta, double sigma2, double nu, bool fixnu);
RcppExport SEXP _ARCensReg_Gradient(SEXP xSEXP, SEXP saemUSEXP, SEXP saemYSEXP, SEXP phiSEXP, SEXP betaSEXP, SEXP sigma2SEXP, SEXP nuSEXP, SEXP fixnuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type saemU(saemUSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type saemY(saemYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type fixnu(fixnuSEXP);
    rcpp_result_gen = Rcpp::wrap(Gradient(x, saemU, saemY, phi, beta, sigma2, nu, fixnu));
    return rcpp_result_gen;
END_RCPP
}
// GradientExp
arma::vec GradientExp(double saemLU, arma::vec saemU, double uy2, arma::vec uyw, arma::mat uw2, arma::mat A, arma::vec B, arma::vec phi, arma::vec beta, double sigma2, double nu, bool fixnu);
RcppExport SEXP _ARCensReg_GradientExp(SEXP saemLUSEXP, SEXP saemUSEXP, SEXP uy2SEXP, SEXP uywSEXP, SEXP uw2SEXP, SEXP ASEXP, SEXP BSEXP, SEXP phiSEXP, SEXP betaSEXP, SEXP sigma2SEXP, SEXP nuSEXP, SEXP fixnuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type saemLU(saemLUSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type saemU(saemUSEXP);
    Rcpp::traits::input_parameter< double >::type uy2(uy2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type uyw(uywSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type uw2(uw2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type fixnu(fixnuSEXP);
    rcpp_result_gen = Rcpp::wrap(GradientExp(saemLU, saemU, uy2, uyw, uw2, A, B, phi, beta, sigma2, nu, fixnu));
    return rcpp_result_gen;
END_RCPP
}
// HessianExp
arma::mat HessianExp(arma::vec saemU, arma::mat saemUZi, arma::vec saemUYi, double uy2, arma::vec uyw, arma::mat uw2, arma::mat A, arma::vec B, arma::vec media, arma::mat X, arma::vec phi, arma::vec beta, double sigma2, double nu, bool fixnu);
RcppExport SEXP _ARCensReg_HessianExp(SEXP saemUSEXP, SEXP saemUZiSEXP, SEXP saemUYiSEXP, SEXP uy2SEXP, SEXP uywSEXP, SEXP uw2SEXP, SEXP ASEXP, SEXP BSEXP, SEXP mediaSEXP, SEXP XSEXP, SEXP phiSEXP, SEXP betaSEXP, SEXP sigma2SEXP, SEXP nuSEXP, SEXP fixnuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type saemU(saemUSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type saemUZi(saemUZiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type saemUYi(saemUYiSEXP);
    Rcpp::traits::input_parameter< double >::type uy2(uy2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type uyw(uywSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type uw2(uw2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type media(mediaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< bool >::type fixnu(fixnuSEXP);
    rcpp_result_gen = Rcpp::wrap(HessianExp(saemU, saemUZi, saemUYi, uy2, uyw, uw2, A, B, media, X, phi, beta, sigma2, nu, fixnu));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ARCensReg_inversa", (DL_FUNC) &_ARCensReg_inversa, 1},
    {"_ARCensReg_ciclobeta", (DL_FUNC) &_ARCensReg_ciclobeta, 5},
    {"_ARCensReg_ComputeMean", (DL_FUNC) &_ARCensReg_ComputeMean, 4},
    {"_ARCensReg_ComputeVar", (DL_FUNC) &_ARCensReg_ComputeVar, 4},
    {"_ARCensReg_Gradient", (DL_FUNC) &_ARCensReg_Gradient, 8},
    {"_ARCensReg_GradientExp", (DL_FUNC) &_ARCensReg_GradientExp, 12},
    {"_ARCensReg_HessianExp", (DL_FUNC) &_ARCensReg_HessianExp, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_ARCensReg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
