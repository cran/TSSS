#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(ngsmthf) (double*, int*, int*, double*, double*, int*, double*, double*, int*, double*, double*, double*, int*, int*, int*, int*);

SEXP NgsmthC(SEXP y, SEXP n, SEXP noisev, SEXP tau2, SEXP bv, SEXP noisew, SEXP sig2, SEXP bw, SEXP initd, SEXP ns, SEXP nfe, SEXP npe, SEXP k1)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8;
    int   *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;

    SEXP ans =  R_NilValue,  trend = R_NilValue, smt = R_NilValue, lkhood = R_NilValue;
    double *xtrend, *xsmt, *xlkhood = NULL;
    int   i, np, nk1;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(noisev);
    d2 = NUMERIC_POINTER(tau2);
    d3 = NUMERIC_POINTER(bv);
    i3 = INTEGER_POINTER(noisew);
    d4 = NUMERIC_POINTER(sig2);
    d5 = NUMERIC_POINTER(bw);
    i4 = INTEGER_POINTER(initd);
    i5 = INTEGER_POINTER(ns);
    i6 = INTEGER_POINTER(nfe);
    i7 = INTEGER_POINTER(npe);
    i8 = INTEGER_POINTER(k1);

    np = *i7;
    nk1 = *i8 + 1;
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, trend = allocVector(REALSXP, np*7));
    SET_VECTOR_ELT(ans, 1, smt = allocVector(REALSXP, nk1*np));
    SET_VECTOR_ELT(ans, 2, lkhood = allocVector(REALSXP, 1)); 

    d6 = NUMERIC_POINTER(trend);
    d7 = NUMERIC_POINTER(smt);
    d8 = NUMERIC_POINTER(lkhood);

    F77_CALL(ngsmthf) (d1,i1,i2,d2,d3,i3,d4,d5,i4,d6,d7,d8,i5,i6,i7,i8);

    xtrend = REAL(trend);
    xsmt = REAL(smt);
    xlkhood = REAL(lkhood);

    for(i=0; i<np*7; i++) xtrend[i] = d6[i];
    for(i=0; i<np*nk1; i++) xsmt[i] = d7[i];
    *xlkhood = *d8;

    UNPROTECT(1);

    return ans;
}
