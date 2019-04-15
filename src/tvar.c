#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(tvarf) (double*, int*, int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*);

SEXP TvarC(SEXP y, SEXP n, SEXP arorder, SEXP torder, SEXP span, SEXP method, SEXP nout, SEXP outlier, SEXP tau20, SEXP delta)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7;

    SEXP ans =  R_NilValue,  tau2 = R_NilValue, sigma2 = R_NilValue, lkhood = R_NilValue, aic = R_NilValue, arcoef = R_NilValue, parcor = R_NilValue;
    double *xtau2, *xsigma2, *xlkhood, *xaic, *xarcoef, *xparcor = NULL;
    int   i, nn, k;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(arorder);
    i3 = INTEGER_POINTER(torder);
    i4 = INTEGER_POINTER(span);
    i5 = INTEGER_POINTER(method);
    i6 = INTEGER_POINTER(nout);
    i7 = INTEGER_POINTER(outlier);
    d2 = NUMERIC_POINTER(tau20);
    d3 = NUMERIC_POINTER(delta);

    nn = (*i1) / (*i4);
    k = (*i2) * nn;
    PROTECT(ans = allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, tau2 = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, sigma2 = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, lkhood = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 3, aic = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 4, arcoef = allocVector(REALSXP, k));
    SET_VECTOR_ELT(ans, 5, parcor = allocVector(REALSXP, k)); 

    d4 = NUMERIC_POINTER(tau2);
    d5 = NUMERIC_POINTER(sigma2);
    d6 = NUMERIC_POINTER(lkhood);
    d7 = NUMERIC_POINTER(aic);
    d8 = NUMERIC_POINTER(arcoef);
    d9 = NUMERIC_POINTER(parcor);

    F77_CALL(tvarf) (d1,i1,i2,i3,i4,i5,i6,i7,d2,d3,d4,d5,d6,d7,d8,d9);

    xtau2 = REAL(tau2);
    xsigma2 = REAL(sigma2);
    xlkhood = REAL(lkhood);
    xaic = REAL(aic);
    xarcoef = REAL(arcoef);
    xparcor = REAL(parcor);

    *xtau2 = *d4;
    *xsigma2 = *d5;
    *xlkhood = *d6;
    *xaic = *d7;
    for(i=0; i<k; i++) xarcoef[i] = d8[i];
    for(i=0; i<k; i++) xparcor[i] = d9[i];

    UNPROTECT(1);

    return ans;
}
