#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(lsar2f)(double*, int*, int*, int*, int*, int*, int*, double*, double*, int*);

SEXP Lsar2C(SEXP y, SEXP n, SEXP k, SEXP n0, SEXP n1, SEXP n2, SEXP ne)
{
    double *d1,*d2,*d3;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7;

    SEXP ans = R_NilValue, aics = R_NilValue, aicmin = R_NilValue, nmin = R_NilValue;
    double *xaics, *xaicmin = NULL;
    int  *xnmin = NULL;
    int m, i;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(k);
    i3 = INTEGER_POINTER(n0);
    i4 = INTEGER_POINTER(n1);
    i5 = INTEGER_POINTER(n2);
    i6 = INTEGER_POINTER(ne);

    m = *i5 - *i4;
    PROTECT(ans = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ans, 0, aics = allocVector(REALSXP, m));
    SET_VECTOR_ELT(ans, 1, aicmin = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 2, nmin = allocVector(INTSXP, 1));

    d2 = NUMERIC_POINTER(aics);
    d3 = NUMERIC_POINTER(aicmin);
    i7 = INTEGER_POINTER(nmin);

    F77_CALL(lsar2f) (d1,i1,i2,i3,i4,i5,i6,d2,d3,i7);

    xaics = REAL(aics);
    xaicmin = REAL(aicmin);
    xnmin = INTEGER(nmin);

    for(i=0; i<m; i++) xaics[i] = d2[i];
    *xaicmin = *d3;
    *xnmin = *i7;

    UNPROTECT(1);

    return ans;
}
