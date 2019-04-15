#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(polregf) (double*,int*,int*,double*,double*,double*,int*,double*);

SEXP PolregC(SEXP y, SEXP n, SEXP k)
{
    double *d1,*d2,*d3,*d4,*d5;
    int *i1,*i2,*i3;

    SEXP ans =  R_NilValue, a = R_NilValue, sig2 = R_NilValue, aic = R_NilValue, imin = R_NilValue, data = R_NilValue;

    double *xa, *xsig2, *xaic, *xdata = NULL;
    int    *ximin= NULL;
    int    i, nn, k1, k2;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(k);

    nn = *i1;
    k1 = *i2 + 1;
    k2 = (*i2) * (*i2);
    PROTECT(ans = allocVector(VECSXP, 5));
    SET_VECTOR_ELT(ans, 0, a = allocVector(REALSXP, k2)); 
    SET_VECTOR_ELT(ans, 1, sig2 = allocVector(REALSXP, k1));
    SET_VECTOR_ELT(ans, 2, aic = allocVector(REALSXP, k1)); 
    SET_VECTOR_ELT(ans, 3, imin = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 4, data = allocVector(REALSXP, nn)); 

    d2 = NUMERIC_POINTER(a);
    d3 = NUMERIC_POINTER(sig2);
    d4 = NUMERIC_POINTER(aic);
    i3 = INTEGER_POINTER(imin);
    d5 = NUMERIC_POINTER(data);

    F77_CALL(polregf) (d1,i1,i2,d2,d3,d4,i3,d5);

    xa = REAL(a);
    xsig2 = REAL(sig2);
    xaic = REAL(aic);
    ximin = INTEGER(imin);
    xdata = REAL(data);

    for(i=0; i<k2; i++) xa[i] = d2[i];
    for(i=0; i<k1; i++) xsig2[i] = d3[i];
    for(i=0; i<k1; i++) xaic[i] = d4[i];
    *ximin = *i3;
    for(i=0; i<nn; i++) xdata[i] = d5[i];

    UNPROTECT(1);

    return ans;
}
