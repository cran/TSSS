#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(crscorf)(double*, int*, int*, int*, double*, double*, double*, double*, double*);

SEXP CrscorC(SEXP y, SEXP n, SEXP l, SEXP lag, SEXP outmin, SEXP outmax)
{
    double *d1, *d2, *d3, *d4, *d5, *d6;
    int *i1, *i2, *i3;

    SEXP ans =  R_NilValue, cov = R_NilValue, cor = R_NilValue, mean = R_NilValue;
    double *xcov, *xcor, *xmean = NULL;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(l);
    i3 = INTEGER_POINTER(lag);
    d2 = NUMERIC_POINTER(outmin);
    d3 = NUMERIC_POINTER(outmax);

    int id = *i2;
    int l2 = id *id;
    int lag1 = *i3+1;
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, cov = allocVector(REALSXP, lag1*l2));
    SET_VECTOR_ELT(ans, 1, cor = allocVector(REALSXP, lag1*l2));
    SET_VECTOR_ELT(ans, 2, mean = allocVector(REALSXP, id)); 

    d4 = NUMERIC_POINTER(cov);
    d5 = NUMERIC_POINTER(cor);
    d6 = NUMERIC_POINTER(mean);

    F77_CALL(crscorf)(d1,i1,i2,i3,d2,d3,d4,d5,d6);

    xcov = REAL(cov);
    xcor = REAL(cor);
    xmean = REAL(mean);

    for(int i=0; i<lag1*l2; i++) xcov[i] = d4[i];
    for(int i=0; i<lag1*l2; i++) xcor[i] = d5[i];
    for(int i=0; i<id; i++) xmean[i] = d6[i];

    UNPROTECT(1);

    return ans;
}
