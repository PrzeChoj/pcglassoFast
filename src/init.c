#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _pcglassoFast_Qinv_down_cpp(void *, void *);
extern SEXP _pcglassoFast_Qinv_down_subet(void *, void *, void *);
extern SEXP _pcglassoFast_Qinv_up_cpp(void *, void *, void *, void *);
extern SEXP _pcglassoFast_updateBeta(void *, void *, void *, void *, void *, void *);
extern SEXP _pcglassoFast_updateLoopCpp(void *, void *, void *, void *, void *, void *, void *);

/* .Fortran calls */
extern void F77_NAME(roptim)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_pcglassoFast_Qinv_down_cpp",   (DL_FUNC) &_pcglassoFast_Qinv_down_cpp,   2},
    {"_pcglassoFast_Qinv_down_subet", (DL_FUNC) &_pcglassoFast_Qinv_down_subet, 3},
    {"_pcglassoFast_Qinv_up_cpp",     (DL_FUNC) &_pcglassoFast_Qinv_up_cpp,     4},
    {"_pcglassoFast_updateBeta",      (DL_FUNC) &_pcglassoFast_updateBeta,      6},
    {"_pcglassoFast_updateLoopCpp",   (DL_FUNC) &_pcglassoFast_updateLoopCpp,   7},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"roptim", (DL_FUNC) &F77_NAME(roptim), 13},
    {NULL, NULL, 0}
};

void R_init_pcglassoFast(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
