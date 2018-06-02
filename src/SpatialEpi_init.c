#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _SpatialEpi_besag_newell_internal(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SpatialEpi_binomialLogLkhd(SEXP, SEXP, SEXP, SEXP);
extern SEXP _SpatialEpi_check_overlap(SEXP, SEXP);
extern SEXP _SpatialEpi_clean_moves_matrix(SEXP, SEXP, SEXP);
extern SEXP _SpatialEpi_coeff(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SpatialEpi_computeAllLogLkhd(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SpatialEpi_kulldorffMC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SpatialEpi_ldmultinom(SEXP, SEXP);
extern SEXP _SpatialEpi_ldnbinom(SEXP, SEXP, SEXP, SEXP);
extern SEXP _SpatialEpi_MCMC_simulation(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SpatialEpi_normalize(SEXP);
extern SEXP _SpatialEpi_NumericVectorEquality(SEXP, SEXP);
extern SEXP _SpatialEpi_poissonLogLkhd(SEXP, SEXP, SEXP, SEXP);
extern SEXP _SpatialEpi_ProbSampleReplace(SEXP);
extern SEXP _SpatialEpi_return_birth_moves(SEXP, SEXP);
extern SEXP _SpatialEpi_return_death_moves(SEXP);
extern SEXP _SpatialEpi_return_local_moves(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_SpatialEpi_besag_newell_internal", (DL_FUNC) &_SpatialEpi_besag_newell_internal, 5},
  {"_SpatialEpi_binomialLogLkhd",       (DL_FUNC) &_SpatialEpi_binomialLogLkhd,       4},
  {"_SpatialEpi_check_overlap",         (DL_FUNC) &_SpatialEpi_check_overlap,         2},
  {"_SpatialEpi_clean_moves_matrix",    (DL_FUNC) &_SpatialEpi_clean_moves_matrix,    3},
  {"_SpatialEpi_coeff",                 (DL_FUNC) &_SpatialEpi_coeff,                 5},
  {"_SpatialEpi_computeAllLogLkhd",     (DL_FUNC) &_SpatialEpi_computeAllLogLkhd,     5},
  {"_SpatialEpi_kulldorffMC",           (DL_FUNC) &_SpatialEpi_kulldorffMC,           5},
  {"_SpatialEpi_ldmultinom",            (DL_FUNC) &_SpatialEpi_ldmultinom,            2},
  {"_SpatialEpi_ldnbinom",              (DL_FUNC) &_SpatialEpi_ldnbinom,              4},
  {"_SpatialEpi_MCMC_simulation",       (DL_FUNC) &_SpatialEpi_MCMC_simulation,       9},
  {"_SpatialEpi_normalize",             (DL_FUNC) &_SpatialEpi_normalize,             1},
  {"_SpatialEpi_NumericVectorEquality", (DL_FUNC) &_SpatialEpi_NumericVectorEquality, 2},
  {"_SpatialEpi_poissonLogLkhd",        (DL_FUNC) &_SpatialEpi_poissonLogLkhd,        4},
  {"_SpatialEpi_ProbSampleReplace",     (DL_FUNC) &_SpatialEpi_ProbSampleReplace,     1},
  {"_SpatialEpi_return_birth_moves",    (DL_FUNC) &_SpatialEpi_return_birth_moves,    2},
  {"_SpatialEpi_return_death_moves",    (DL_FUNC) &_SpatialEpi_return_death_moves,    1},
  {"_SpatialEpi_return_local_moves",    (DL_FUNC) &_SpatialEpi_return_local_moves,    3},
  {NULL, NULL, 0}
};

void R_init_SpatialEpi(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}