#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(periodf)(double*, int*, int*, int*, double*, double*, double*, double*, int*);

SEXP PeriodC(SEXP y, SEXP n, SEXP np, SEXP iw, SEXP outmin, SEXP outmax)
{
    double *d1, *d2, *d3, *d4, *d5;
    int *i1, *i2, *i3, *i4;
    int np1;

    SEXP ans =  R_NilValue, pe = R_NilValue, spe = R_NilValue, ifg = R_NilValue;
    double *xpe, *xspe = NULL;
    int  *xifg;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(np);
    i3 = INTEGER_POINTER(iw);
    d2 = NUMERIC_POINTER(outmin);
    d3 = NUMERIC_POINTER(outmax);
    np1 = *i2 + 1;


    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, pe = allocVector(REALSXP, np1));
    SET_VECTOR_ELT(ans, 1, spe = allocVector(REALSXP, np1));
    SET_VECTOR_ELT(ans, 2, np = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 3, ifg = allocVector(INTSXP, 1));

    d4 = NUMERIC_POINTER(pe);
    d5 = NUMERIC_POINTER(spe);
    i4 = INTEGER_POINTER(ifg);

    F77_CALL(periodf)(d1,i1,i2,i3,d2,d3,d4,d5,i4);

    xpe = REAL(pe);
    xspe = REAL(spe);
    xifg = INTEGER(ifg);

    for(int i=0; i<np1; i++) xpe[i] = d4[i];
    for(int i=0; i<np1; i++) xspe[i] = d5[i];
    *xifg = *i4;

    UNPROTECT(1);

    return ans;
}
