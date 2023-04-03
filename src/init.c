#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> //
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/



/* .C calls */

extern void GAUSSIAN_PPMX(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                          void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                          void *, void *, void *, void *, void *,void *);

extern void ORDINAL_PPMX(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                         void *, void *, void *, void *, void *,void *, void *, void *);

extern void GAUSSIAN_PPMX_MISSING(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                  void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                  void *, void *, void *, void *, void *,void *, void *, void *, void *, void *);

extern void ORDINAL_PPMX_MISSING(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                 void *, void *);

extern void mcmc_sppm(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                      void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                      void *, void *, void *, void *, void *,void *, void *);


extern void rppmx(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                  void *, void *, void *, void *, void *,void *, void *, void *, void *, void *);


extern void mcmc_curvecluster(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                              void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                              void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                              void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                              void *, void *, void *, void *, void *,void *, void *);

extern void mcmc_bivariate_curvecluster(void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                        void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                        void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                        void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                        void *, void *, void *, void *, void *,void *, void *, void *, void *, void *,
                                        void *, void *, void *, void *, void *,void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"gaussian_ppmx", (DL_FUNC) &GAUSSIAN_PPMX, 26},
    {"ordinal_ppmx", (DL_FUNC) &ORDINAL_PPMX, 28},
    {"gaussian_ppmx_missing", (DL_FUNC) &GAUSSIAN_PPMX_MISSING, 30},
    {"ordinal_ppmx_missing", (DL_FUNC) &ORDINAL_PPMX_MISSING, 32},
    {"mcmc_sppm", (DL_FUNC) &mcmc_sppm, 27},
    {"rppmx", (DL_FUNC) &rppmx, 20},
    {"mcmc_curvecluster", (DL_FUNC) &mcmc_curvecluster, 47},
    {"mcmc_bivariate_curvecluster", (DL_FUNC) &mcmc_bivariate_curvecluster, 58},
    {NULL, NULL, 0}
};

void R_init_modernVA(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
