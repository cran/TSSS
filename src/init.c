#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP ArfitC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP arma(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP armaft(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP arsp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP BoxcoxC(SEXP, SEXP);
extern SEXP CrscorC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP densty(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FftperC(SEXP, SEXP, SEXP);
extern SEXP KlinfoC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Lsar1C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Lsar2C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LsqrC(SEXP, SEXP, SEXP, SEXP);
extern SEXP MarfitC(SEXP, SEXP, SEXP, SEXP);
extern SEXP MarlsqC(SEXP, SEXP, SEXP, SEXP);
extern SEXP MarspcC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NgsimC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP NgsmthC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PeriodC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PolregC(SEXP, SEXP, SEXP);
extern SEXP SeasonC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP SimssmC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP smooth(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TrendC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TvarC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TvspcC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TvvarC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP UnicorC(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"ArfitC",  (DL_FUNC) &ArfitC,   6},
    {"arma",    (DL_FUNC) &arma,     9},
    {"armaft",  (DL_FUNC) &armaft,   8},
    {"arsp",    (DL_FUNC) &arsp,     6},
    {"BoxcoxC", (DL_FUNC) &BoxcoxC,  2},
    {"CrscorC", (DL_FUNC) &CrscorC,  6},
    {"densty",  (DL_FUNC) &densty,   5},
    {"FftperC", (DL_FUNC) &FftperC,  3},
    {"KlinfoC", (DL_FUNC) &KlinfoC,  6},
    {"Lsar1C",  (DL_FUNC) &Lsar1C,   6},
    {"Lsar2C",  (DL_FUNC) &Lsar2C,   7},
    {"LsqrC",   (DL_FUNC) &LsqrC,    4},
    {"MarfitC", (DL_FUNC) &MarfitC,  4},
    {"MarlsqC", (DL_FUNC) &MarlsqC,  4},
    {"MarspcC", (DL_FUNC) &MarspcC,  5},
    {"NgsimC",  (DL_FUNC) &NgsimC,  18},
    {"NgsmthC", (DL_FUNC) &NgsmthC, 13},
    {"PeriodC", (DL_FUNC) &PeriodC,  6},
    {"PolregC", (DL_FUNC) &PolregC,  3},
    {"SeasonC", (DL_FUNC) &SeasonC, 20},
    {"SimssmC", (DL_FUNC) &SimssmC, 14},
    {"smooth",  (DL_FUNC) &smooth,  18},
    {"TrendC",  (DL_FUNC) &TrendC,   6},
    {"TvarC",   (DL_FUNC) &TvarC,   10},
    {"TvspcC",  (DL_FUNC) &TvspcC,   8},
    {"TvvarC",  (DL_FUNC) &TvvarC,   6},
    {"UnicorC", (DL_FUNC) &UnicorC,  5},
    {NULL, NULL, 0}
};

void R_init_TSSS(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
