#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(pfilternf)(double*, int*, int*, int*, double*, double*, double*, double*, int*, double*, double*);

SEXP PfiltnlC(SEXP y, SEXP n, SEXP m, SEXP lag, SEXP sig2, SEXP tau2, SEXP xmin, SEXP xmax, SEXP ix)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7;
    int *i1,*i2,*i3,*i4;

    SEXP ans = R_NilValue, t = R_NilValue, ff = R_NilValue;
    double *xt, *xff = NULL;
    int   i, nn;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(m);
    i3 = INTEGER_POINTER(lag);
    d2 = NUMERIC_POINTER(sig2);
    d3 = NUMERIC_POINTER(tau2);
    d4 = NUMERIC_POINTER(xmin);
    d5 = NUMERIC_POINTER(xmax);
    i4 = INTEGER_POINTER(ix);

    nn = *i1;
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, t = allocVector(REALSXP, 8*nn));
    SET_VECTOR_ELT(ans, 1, ff = allocVector(REALSXP, 1));

    d6 = NUMERIC_POINTER(t);
    d7 = NUMERIC_POINTER(ff);

    F77_CALL(pfilternf) (d1,i1,i2,i3,d2,d3,d4,d5,i4,d6,d7);

    xt = REAL(t);
    xff = REAL(ff);

    for(i=0; i<8*nn; i++) xt[i] = d6[i];
    *xff = *d7;

    UNPROTECT(1);

    return ans;
}

