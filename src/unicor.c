#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(unicorf)(double*, int*, int*, double*, double*, double*, double*);

SEXP UnicorC(SEXP y, SEXP n, SEXP lag, SEXP outmin, SEXP outmax)
{
    double *d1, *d2, *d3, *d4, *d5;
    int *i1, *i2;

    SEXP ans =  R_NilValue, cov = R_NilValue, mean = R_NilValue;
    double *xcov, *xmean = NULL;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(lag);
    d2 = NUMERIC_POINTER(outmin);
    d3 = NUMERIC_POINTER(outmax);

    int lag1 = *i2+1;
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, cov = allocVector(REALSXP, lag1*4));
    SET_VECTOR_ELT(ans, 1, mean = allocVector(REALSXP, 1)); 

    d4 = NUMERIC_POINTER(cov);
    d5 = NUMERIC_POINTER(mean);

    F77_CALL(unicorf)(d1,i1,i2,d2,d3,d4,d5);

    xcov = REAL(cov);
    xmean = REAL(mean);

    for(int i=0; i<lag1*4; i++) xcov[i] = d4[i];
    *xmean = *d5;

    UNPROTECT(1);

    return ans;
}
