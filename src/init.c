#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _dFCaHSMM_ck_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_covmat_c(SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_create_B_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_create_Bii_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_create_Bij_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_create_gradB_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_create_gradBii_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_create_gradBij_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_dmvnrm_arma(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_em_estep_covar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_em_mstep_covar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_Fk_cpp(SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_forward_backward_hsmm_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_get_Ahat(SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_get_emission_distribution(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_get_lambdan(SEXP, SEXP);
extern SEXP _dFCaHSMM_get_locations(SEXP);
extern SEXP _dFCaHSMM_grad_ck_cpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_grad_term2_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_hsmm_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_logsumexp_cpp(SEXP);
extern SEXP _dFCaHSMM_logvec_c(SEXP);
extern SEXP _dFCaHSMM_optim_covar_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_params2vec_covar(SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_pk_cpp(SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_standardize_rows_cpp(SEXP);
extern SEXP _dFCaHSMM_term2_covar_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dFCaHSMM_vec2params_covar(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_dFCaHSMM_ck_cpp",                    (DL_FUNC) &_dFCaHSMM_ck_cpp,                     4},
  {"_dFCaHSMM_covmat_c",                  (DL_FUNC) &_dFCaHSMM_covmat_c,                   3},
  {"_dFCaHSMM_create_B_cpp",              (DL_FUNC) &_dFCaHSMM_create_B_cpp,               5},
  {"_dFCaHSMM_create_Bii_cpp",            (DL_FUNC) &_dFCaHSMM_create_Bii_cpp,             4},
  {"_dFCaHSMM_create_Bij_cpp",            (DL_FUNC) &_dFCaHSMM_create_Bij_cpp,             5},
  {"_dFCaHSMM_create_gradB_cpp",          (DL_FUNC) &_dFCaHSMM_create_gradB_cpp,           4},
  {"_dFCaHSMM_create_gradBii_cpp",        (DL_FUNC) &_dFCaHSMM_create_gradBii_cpp,         4},
  {"_dFCaHSMM_create_gradBij_cpp",        (DL_FUNC) &_dFCaHSMM_create_gradBij_cpp,         4},
  {"_dFCaHSMM_dmvnrm_arma",               (DL_FUNC) &_dFCaHSMM_dmvnrm_arma,                5},
  {"_dFCaHSMM_em_estep_covar",            (DL_FUNC) &_dFCaHSMM_em_estep_covar,            14},
  {"_dFCaHSMM_em_mstep_covar",            (DL_FUNC) &_dFCaHSMM_em_mstep_covar,            14},
  {"_dFCaHSMM_Fk_cpp",                    (DL_FUNC) &_dFCaHSMM_Fk_cpp,                     3},
  {"_dFCaHSMM_forward_backward_hsmm_cpp", (DL_FUNC) &_dFCaHSMM_forward_backward_hsmm_cpp,  7},
  {"_dFCaHSMM_get_Ahat",                  (DL_FUNC) &_dFCaHSMM_get_Ahat,                   3},
  {"_dFCaHSMM_get_emission_distribution", (DL_FUNC) &_dFCaHSMM_get_emission_distribution,  7},
  {"_dFCaHSMM_get_lambdan",               (DL_FUNC) &_dFCaHSMM_get_lambdan,                2},
  {"_dFCaHSMM_get_locations",             (DL_FUNC) &_dFCaHSMM_get_locations,              1},
  {"_dFCaHSMM_grad_ck_cpp",               (DL_FUNC) &_dFCaHSMM_grad_ck_cpp,                4},
  {"_dFCaHSMM_grad_term2_cpp",            (DL_FUNC) &_dFCaHSMM_grad_term2_cpp,             5},
  {"_dFCaHSMM_hsmm_cpp",                  (DL_FUNC) &_dFCaHSMM_hsmm_cpp,                  11},
  {"_dFCaHSMM_logsumexp_cpp",             (DL_FUNC) &_dFCaHSMM_logsumexp_cpp,              1},
  {"_dFCaHSMM_logvec_c",                  (DL_FUNC) &_dFCaHSMM_logvec_c,                   1},
  {"_dFCaHSMM_optim_covar_rcpp",          (DL_FUNC) &_dFCaHSMM_optim_covar_rcpp,           5},
  {"_dFCaHSMM_params2vec_covar",          (DL_FUNC) &_dFCaHSMM_params2vec_covar,           3},
  {"_dFCaHSMM_pk_cpp",                    (DL_FUNC) &_dFCaHSMM_pk_cpp,                     3},
  {"_dFCaHSMM_standardize_rows_cpp",      (DL_FUNC) &_dFCaHSMM_standardize_rows_cpp,       1},
  {"_dFCaHSMM_term2_covar_cpp",           (DL_FUNC) &_dFCaHSMM_term2_covar_cpp,            6},
  {"_dFCaHSMM_vec2params_covar",          (DL_FUNC) &_dFCaHSMM_vec2params_covar,           4},
  {NULL, NULL, 0}
};

void R_init_dFCaHSMM(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
