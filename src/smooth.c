#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(smoothf) (double*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, int*, int*, double*, double*, int*, int*, int*, double*, double*, double*, double*);

SEXP smooth(SEXP yy, SEXP n, SEXP m, SEXP k, SEXP ff, SEXP gg, SEXP hh, SEXP qq, SEXP rr, SEXP x0, SEXP v0, SEXP fend, SEXP pend, SEXP omin, SEXP omax, SEXP nmiss, SEXP startp, SEXP np)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;

    SEXP ans = R_NilValue, xss = R_NilValue, vss = R_NilValue, lkhood = R_NilValue, aic = R_NilValue;
    double *xxss, *xvss, *xlkhood, *xaic = NULL;
    int   i, mm,npe;

    d1 = NUMERIC_POINTER(yy);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(m);
    i3 = INTEGER_POINTER(k);
    d2 = NUMERIC_POINTER(ff);
    d3 = NUMERIC_POINTER(gg);
    d4 = NUMERIC_POINTER(hh);
    d5 = NUMERIC_POINTER(qq);
    d6 = NUMERIC_POINTER(rr);
    d7 = NUMERIC_POINTER(x0);
    d8 = NUMERIC_POINTER(v0);
    i4 = INTEGER_POINTER(fend);
    i5 = INTEGER_POINTER(pend);
    d9 = NUMERIC_POINTER(omin);
    d10 = NUMERIC_POINTER(omax);
    i6 = INTEGER_POINTER(nmiss);
    i7 = INTEGER_POINTER(startp);
    i8 = INTEGER_POINTER(np);

    mm = *i2;
    npe = *i5;
    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, xss = allocVector(REALSXP, mm*npe));
    SET_VECTOR_ELT(ans, 1, vss = allocVector(REALSXP, mm*mm*npe));
    SET_VECTOR_ELT(ans, 2, lkhood = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 3, aic = allocVector(REALSXP, 1));

    d11 = NUMERIC_POINTER(xss);
    d12 = NUMERIC_POINTER(vss);
    d13 = NUMERIC_POINTER(lkhood);
    d14 = NUMERIC_POINTER(aic);

    F77_CALL(smoothf) (d1,i1,i2,i3,d2,d3,d4,d5,d6,d7,d8,i4,i5,d9,d10,i6,i7,i8,d11,d12,d13,d14);

    xxss = REAL(xss);
    xvss = REAL(vss);
    xlkhood = REAL(lkhood);
    xaic = REAL(aic);

    for(i=0; i<mm*npe; i++) xxss[i] = d11[i];
    for(i=0; i<mm*mm*npe; i++) xvss[i] = d12[i];
    *xlkhood = *d13;
    *xaic = *d14;

    UNPROTECT(1);

    return ans;
}

