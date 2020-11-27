#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> //
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/



/* .C calls */
extern void mcmc_missing(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                         void *, void *, void *);


extern void ordinal_missing_ppmx(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,void *, void *, void *, void *);


extern void ordinal_ppmx(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *);

extern void mcmcppmx(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                     void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                     void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                     void *, void *, void *, void *, void *,void *, void *, void *, void *);

extern void mcmc_sppm(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                      void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                      void *, void *, void *, void *, void *,void *, void *);


extern void rppmx(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                  void *, void *, void *, void *, void *,void *, void *, void *, void *, void *);


static const R_CMethodDef CEntries[] = {
    {"mcmc_missing", (DL_FUNC) &mcmc_missing, 43},
    {"ordinal_missing_ppmx", (DL_FUNC) &ordinal_missing_ppmx, 49},
    {"ordinal_ppmx", (DL_FUNC) &ordinal_ppmx, 45},
    {"mcmcppmx", (DL_FUNC) &mcmcppmx, 39},
    {"mcmc_sppm", (DL_FUNC) &mcmc_sppm, 27},
    {"rppmx", (DL_FUNC) &rppmx, 20},
    {NULL, NULL, 0}
};

void R_init_modernVA(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
