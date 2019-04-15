#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(simssmf)  (int*,int*,int*,int*,int*,int*,int*,double*,int*,double*,double*,double*,double*,double*,double*);

SEXP SimssmC(SEXP m1, SEXP m2, SEXP m3, SEXP m, SEXP k, SEXP n, SEXP ini, SEXP sig2, SEXP period, SEXP tau1, SEXP tau2, SEXP tau3, SEXP arcoef, SEXP x)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;

    SEXP ans = R_NilValue, y = R_NilValue;
    double *yy= NULL;
    int   i, nn;

    i1 = INTEGER_POINTER(m1);
    i2 = INTEGER_POINTER(m2);
    i3 = INTEGER_POINTER(m3);
    i4 = INTEGER_POINTER(m);
    i5 = INTEGER_POINTER(k);
    i6 = INTEGER_POINTER(n);
    i7 = INTEGER_POINTER(ini);
    d1 = NUMERIC_POINTER(sig2);
    i8 = INTEGER_POINTER(period);
    d2 = NUMERIC_POINTER(tau1);
    d3 = NUMERIC_POINTER(tau2);
    d4 = NUMERIC_POINTER(tau3);
    d5 = NUMERIC_POINTER(arcoef);
    d6 = NUMERIC_POINTER(x);

    nn = *i6;
    PROTECT(ans = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(ans, 0, y = allocVector(REALSXP, nn));

    d7 = NUMERIC_POINTER(y);

    F77_CALL(simssmf) (i1,i2,i3,i4,i5,i6,i7,d1,i8,d2,d3,d4,d5,d6,d7);

    yy = REAL(y);

    for(i=0; i<nn; i++) yy[i] = d7[i];
    UNPROTECT(1);

    return ans;
}

