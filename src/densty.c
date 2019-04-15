#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(denstyf) ( int*, double*, double*, double*, int*, double* );

SEXP densty(SEXP model, SEXP param, SEXP xmin, SEXP xmax, SEXP k)
{
    double *d1,*d2,*d3,*d4;
    int   *i1,*i2;

    SEXP ans =  R_NilValue,  f = R_NilValue;
    double *xf = NULL;
    int   i, nk;

    i1 = INTEGER_POINTER(model);
    d1 = NUMERIC_POINTER(param);
    d2 = NUMERIC_POINTER(xmin);
    d3 = NUMERIC_POINTER(xmax);
    i2 = INTEGER_POINTER(k);

    nk = *i2;
    PROTECT(ans = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(ans, 0, f = allocVector(REALSXP, nk));

    d4 = NUMERIC_POINTER(f);

    F77_CALL(denstyf) (i1,d1,d2,d3,i2,d4);

    xf = REAL(f);

    for(i=0; i<nk; i++) xf[i] = d4[i];

    UNPROTECT(1);

    return ans;
}
