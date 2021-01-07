#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(trendf) (double*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*);

SEXP TrendC(SEXP y, SEXP n, SEXP m, SEXP iopt, SEXP tau20, SEXP delta)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int *i1,*i2,*i3;

    SEXP ans =  R_NilValue, tau2 = R_NilValue, sig2 = R_NilValue, lkhood = R_NilValue, aic = R_NilValue, xss =  R_NilValue, rs = R_NilValue;

    double *xtau2, *xsig2, *xlkhood, *xaic, *xxss, *xrs = NULL;
    int   i, mm, nn;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(m);
    i3 = INTEGER_POINTER(iopt);
    d2 = NUMERIC_POINTER(tau20);
    d3 = NUMERIC_POINTER(delta);

    nn = *i1;
    mm = *i2;
    PROTECT(ans = allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, tau2 = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 1, sig2 = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, lkhood = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 3, aic = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 4, xss = allocVector(REALSXP, mm*nn));
    SET_VECTOR_ELT(ans, 5, rs = allocVector(REALSXP, nn)); 

    d4 = NUMERIC_POINTER(tau2);
    d5 = NUMERIC_POINTER(sig2);
    d6 = NUMERIC_POINTER(lkhood);
    d7 = NUMERIC_POINTER(aic);
    d8 = NUMERIC_POINTER(xss);
    d9 = NUMERIC_POINTER(rs);

    F77_CALL(trendf) (d1,i1,i2,i3,d2,d3,d4,d5,d6,d7,d8,d9);

    xtau2 = REAL(tau2);
    xsig2 = REAL(sig2);
    xlkhood = REAL(lkhood);
    xaic = REAL(aic);
    xxss = REAL(xss);
    xrs = REAL(rs);

    *xtau2 = *d4;
    *xsig2 = *d5;
    *xlkhood = *d6;
    *xaic = *d7;
    for(i=0; i<mm*nn; i++) xxss[i] = d8[i];
    for(i=0; i<nn; i++) xrs[i] = d9[i];

    UNPROTECT(1);

    return ans;
}

