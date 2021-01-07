#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(armaft2f)(double*, int*, int*, int*, int*, double*, double*, double*, double*, double*, int*);

SEXP armaft2(SEXP y, SEXP n, SEXP mmax, SEXP lmax, SEXP mlmax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6;
    int *i1,*i2,*i3,*i4,*i5;

    SEXP ans = R_NilValue, sig2 = R_NilValue, flk = R_NilValue, aic = R_NilValue, ar = R_NilValue,  ma = R_NilValue, ier = R_NilValue;
    double *xsig2, *xflk, *xaic, *xar, *xma = NULL;
    int *xier = NULL;
    int mm, ll, mm1, ll1;
    int i;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(mmax);
    i3 = INTEGER_POINTER(lmax);
    i4 = INTEGER_POINTER(mlmax);

    mm = *i2;
    ll = *i3;
    mm1 = mm + 1;
    ll1 = ll + 1;

    PROTECT(ans = allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, sig2 = allocVector(REALSXP, mm1*ll1));
    SET_VECTOR_ELT(ans, 1, flk = allocVector(REALSXP, mm1*ll1));
    SET_VECTOR_ELT(ans, 2, aic = allocVector(REALSXP, mm1*ll1));
    SET_VECTOR_ELT(ans, 3, ar = allocVector(REALSXP, mm1*ll1*mm));
    SET_VECTOR_ELT(ans, 4, ma = allocVector(REALSXP, mm1*ll1*ll));
    SET_VECTOR_ELT(ans, 5, ier = allocVector(INTSXP, 1));

    d2 = NUMERIC_POINTER(sig2);
    d3 = NUMERIC_POINTER(flk);
    d4 = NUMERIC_POINTER(aic);
    d5 = NUMERIC_POINTER(ar);
    d6 = NUMERIC_POINTER(ma);
    i5 = INTEGER_POINTER(ier);

    F77_CALL(armaft2f) (d1,i1,i2,i3,i4,d2,d3,d4,d5,d6,i5);

    xsig2 = REAL(sig2);
    xflk = REAL(flk);
    xaic = REAL(aic);
    xar = REAL(ar);
    xma = REAL(ma);
    xier = INTEGER(ier);

    for(i=0; i<mm1*ll1; i++) xsig2[i] = d2[i];
    for(i=0; i<mm1*ll1; i++) xflk[i] = d3[i];
    for(i=0; i<mm1*ll1; i++) xaic[i] = d4[i];
    for(i=0; i<mm1*ll1*mm; i++) xar[i] = d5[i];
    for(i=0; i<mm1*ll1*ll; i++) xma[i] = d6[i];
    *xier = *i5;

    UNPROTECT(1);

    return ans;
}
