#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(lsqrf)(double*, int*, int*, int*, double*, double*, int*, double*, double*);

SEXP LsqrC(SEXP y, SEXP n, SEXP k, SEXP mj)
{
    double *d1,*d2,*d3, *d4, *d5;
    int *i1,*i2,*i3,*i4;

    SEXP ans = R_NilValue, aic = R_NilValue, sig2 = R_NilValue, imin = R_NilValue, reg = R_NilValue, data = R_NilValue;
    double *xaic, *xsig2, *xreg, *xdata = NULL;
    int  *ximin = NULL;
    int  k1, kk, nn, i;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(k);
    i3 = INTEGER_POINTER(mj);

    kk = (*i2) * (*i2);
    k1 = *i2 +1;
    nn = *i1;
    PROTECT(ans = allocVector(VECSXP, 5));
    SET_VECTOR_ELT(ans, 0, aic = allocVector(REALSXP, k1));
    SET_VECTOR_ELT(ans, 1, sig2 = allocVector(REALSXP, k1));
    SET_VECTOR_ELT(ans, 2, imin = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 3, reg = allocVector(REALSXP, kk));
    SET_VECTOR_ELT(ans, 4, data = allocVector(REALSXP, nn));

    d2 = NUMERIC_POINTER(aic);
    d3 = NUMERIC_POINTER(sig2);
    i4 = INTEGER_POINTER(imin);
    d4 = NUMERIC_POINTER(reg);
    d5 = NUMERIC_POINTER(data);

    F77_CALL(lsqrf) (d1,i1,i2,i3,d2,d3,i4,d4,d5);

    xaic = REAL(aic);
    xsig2 = REAL(sig2);
    ximin = INTEGER(imin);
    xreg = REAL(reg);
    xdata = REAL(data);

    for(i=0; i<k1; i++) xaic[i] = d2[i];
    for(i=0; i<k1; i++) xsig2[i] = d3[i];
    *ximin = *i4;
    for(i=0; i<kk; i++) xreg[i] = d4[i];
    for(i=0; i<nn; i++) xdata[i] = d5[i];

    UNPROTECT(1);

    return ans;
}
