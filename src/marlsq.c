#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(marlsqf)(double*, int*, int*, int*, double*, double*, int*, double*);

SEXP MarlsqC(SEXP y, SEXP n, SEXP l, SEXP lag)
{
    double *d1, *d2, *d3, *d4;
    int      *i1, *i2, *i3, *i4;

    SEXP ans = R_NilValue, arcoef = R_NilValue, v = R_NilValue, lmax = R_NilValue, aic = R_NilValue;
    double *xarcoef, *xv, *xaic = NULL;
    int  *xlmax = NULL;
    int  i, ll, lg;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(l);
    i3 = INTEGER_POINTER(lag);

    ll = (*i2) * (*i2);
    lg = *i3;
    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, arcoef = allocVector(REALSXP, ll*lg));
    SET_VECTOR_ELT(ans, 1, v = allocVector(REALSXP, ll));
    SET_VECTOR_ELT(ans, 2, lmax = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 3, aic = allocVector(REALSXP, 1));

    d2 = NUMERIC_POINTER(arcoef);
    d3 = NUMERIC_POINTER(v);
    i4 = INTEGER_POINTER(lmax);
    d4 = NUMERIC_POINTER(aic);

    F77_CALL(marlsqf) (d1,i1,i2,i3,d2,d3,i4,d4);

    xarcoef = REAL(arcoef);
    xv = REAL(v);
    xlmax = INTEGER(lmax);
    xaic = REAL(aic);

    for(i=0; i<ll*lg; i++) xarcoef[i] = d2[i];
    for(i=0; i<ll; i++) xv[i] = d3[i];
    *xlmax = *i4;
    *xaic = *d4;

    UNPROTECT(1);

    return ans;
}
