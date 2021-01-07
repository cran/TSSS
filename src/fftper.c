#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(fftperf)(double*, int*, int*, double*, double*, int*, int*);

SEXP FftperC(SEXP y, SEXP n, SEXP iw)
{
    double *d1, *d2, *d3;
    int *i1, *i2, *i3, *i4;

    SEXP ans =  R_NilValue, pe = R_NilValue, spe = R_NilValue, np = R_NilValue, ifg = R_NilValue;
    double *xpe, *xspe = NULL;
    int *xnp, *xifg = NULL;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(iw);

    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, pe = allocVector(REALSXP, 1025));
    SET_VECTOR_ELT(ans, 1, spe = allocVector(REALSXP, 1025));
    SET_VECTOR_ELT(ans, 2, np = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 3, ifg = allocVector(INTSXP, 1));

    d2 = NUMERIC_POINTER(pe);
    d3 = NUMERIC_POINTER(spe);
    i3 = INTEGER_POINTER(np);
    i4 = INTEGER_POINTER(ifg);

    F77_CALL(fftperf)(d1,i1,i2,d2,d3,i3,i4);

    xpe = REAL(pe);
    xspe = REAL(spe);
    xnp = INTEGER(np);
    xifg = INTEGER(ifg);

    for(int i=0; i<1025; i++) xpe[i] = d2[i];
    for(int i=0; i<1025; i++) xspe[i] = d3[i];
    *xnp = *i3;
    *xifg = *i4;

    UNPROTECT(1);

    return ans;
}
