#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(tvvarf) (double*, int*, int*, double*, int*, double*, double*, double*, double*, int*, double*, double*, double*, double*, double*, double*);

SEXP TvvarC(SEXP y, SEXP n, SEXP m, SEXP tau20, SEXP iopt, SEXP delta)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12;
    int *i1,*i2,*i3,*i4;

    SEXP ans =  R_NilValue,  tvvar = R_NilValue, normdat = R_NilValue, y1 = R_NilValue,  n1 = R_NilValue, trend = R_NilValue, noise = R_NilValue, taumax = R_NilValue, sig2m = R_NilValue, ffmax = R_NilValue, aic = R_NilValue;
    double *xtvvar, *xnrdat, *xy1, *xtrend, *xnoise, *xtaumax, *xsig2m, *xffmax, *xaic = NULL;
    int    *xn1;
    int   i, nn, n2;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(m);
    d2 = NUMERIC_POINTER(tau20);
    i3 = INTEGER_POINTER(iopt);
    d3 = NUMERIC_POINTER(delta);

    nn = *i1;
    n2 = nn/2;
    PROTECT(ans = allocVector(VECSXP, 10));
    SET_VECTOR_ELT(ans, 0, tvvar = allocVector(REALSXP, n2));
    SET_VECTOR_ELT(ans, 1, normdat = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 2, y1 = allocVector(REALSXP, n2)); 
    SET_VECTOR_ELT(ans, 3, n1 = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 4, trend = allocVector(REALSXP, n2*3));
    SET_VECTOR_ELT(ans, 5, noise = allocVector(REALSXP, n2)); 
    SET_VECTOR_ELT(ans, 6, taumax = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 7, sig2m = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 8, ffmax = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 9, aic = allocVector(REALSXP, 1)); 

    d4 = NUMERIC_POINTER(tvvar);
    d5 = NUMERIC_POINTER(normdat);
    d6 = NUMERIC_POINTER(y1);
    i4 = INTEGER_POINTER(n1);
    d7 = NUMERIC_POINTER(trend);
    d8 = NUMERIC_POINTER(noise);
    d9 = NUMERIC_POINTER(taumax);
    d10 = NUMERIC_POINTER(sig2m);
    d11 = NUMERIC_POINTER(ffmax);
    d12 = NUMERIC_POINTER(aic);

    F77_CALL(tvvarf) (d1,i1,i2,d2,i3,d3,d4,d5,d6,i4,d7,d8,d9,d10,d11,d12);

    xtvvar = REAL(tvvar);
    xnrdat = REAL(normdat);
    xy1 = REAL(y1);
    xn1 = INTEGER(n1);
    xtrend = REAL(trend);
    xnoise = REAL(noise);
    xtaumax = REAL(taumax);
    xsig2m = REAL(sig2m);
    xffmax = REAL(ffmax);
    xaic = REAL(aic);

    for(i=0; i<n2; i++) xtvvar[i] = d4[i];
    for(i=0; i<nn; i++) xnrdat[i] = d5[i];
    for(i=0; i<n2; i++) xy1[i] = d6[i];
    *xn1 = *i4;
    for(i=0; i<n2*3; i++) xtrend[i] = d7[i];
    for(i=0; i<n2; i++) xnoise[i] = d8[i];
    *xtaumax = *d9;
    *xsig2m = *d10;
    *xffmax = *d11;
    *xaic = *d12;

    UNPROTECT(1);

    return ans;
}

