#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(armaftf)(double*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*,  double*, double*, int*);

SEXP armaft(SEXP y, SEXP n, SEXP m, SEXP l, SEXP mlmax, SEXP iparam, SEXP ar0, SEXP ma0)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8;
    int *i1,*i2,*i3,*i4,*i5,*i6;

    SEXP ans = R_NilValue, sig2 = R_NilValue, flk = R_NilValue, aic = R_NilValue, ar = R_NilValue,  ma = R_NilValue, ier = R_NilValue;
    double *xsig2, *xflk, *xaic, *xar, *xma = NULL;
    int *xier = NULL;
    int mm, ll;
    int i;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(m);
    i3 = INTEGER_POINTER(l);
    i4 = INTEGER_POINTER(mlmax);
    i5 = INTEGER_POINTER(iparam);
    d2 = NUMERIC_POINTER(ar0);
    d3 = NUMERIC_POINTER(ma0);

    mm = *i2;
    ll = *i3;

    PROTECT(ans = allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, sig2 = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, flk = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, aic = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 3, ar = allocVector(REALSXP, mm));
    SET_VECTOR_ELT(ans, 4, ma = allocVector(REALSXP, ll));
    SET_VECTOR_ELT(ans, 5, ier = allocVector(INTSXP, 1));

    d4 = NUMERIC_POINTER(sig2);
    d5 = NUMERIC_POINTER(flk);
    d6 = NUMERIC_POINTER(aic);
    d7 = NUMERIC_POINTER(ar);
    d8 = NUMERIC_POINTER(ma);
    i6 = INTEGER_POINTER(ier);

    F77_CALL(armaftf) (d1,i1,i2,i3,i4,i5,d2,d3,d4,d5,d6,d7,d8,i6);

    xsig2 = REAL(sig2);
    xflk = REAL(flk);
    xaic = REAL(aic);
    xar = REAL(ar);
    xma = REAL(ma);
    xier = INTEGER(ier);

    *xsig2 = *d4;
    *xflk = *d5;
    *xaic = *d6;
    for(i=0; i<mm; i++) xar[i] = d7[i];
    for(i=0; i<ll; i++) xma[i] = d8[i];
    *xier = *i6;

    UNPROTECT(1);

    return ans;
}
