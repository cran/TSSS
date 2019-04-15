#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(seasonf) (double*, int*, int*, int*, int*, int*, int*, int*, int*, double*, int*, int*, int*, double*, int*, int*,  double*, double*, int*, int*, double*, double*, double*, double*, double*, double*, int*, int*);

SEXP SeasonC(SEXP y, SEXP n, SEXP m1, SEXP m2, SEXP m3, SEXP m4, SEXP iper, SEXP jyear, SEXP month, SEXP tau2, SEXP ns, SEXP nfe, SEXP npe, SEXP ar, SEXP logt, SEXP iopt, SEXP omin, SEXP omax, SEXP nmax, SEXP mj)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11;
    int *i1,*i2,*i3, *i4, *i5, *i6, *i7, *i8, *i9, *i10, *i11, *i12, *i13, *i14, *i15, *i16, *i17;

    SEXP ans =  R_NilValue, ff = R_NilValue, over = R_NilValue, aic = R_NilValue, xss =  R_NilValue, vss = R_NilValue, deff = R_NilValue, ier1 = R_NilValue, ier2 = R_NilValue;

    double *xff, *xover, *xaic, *xxss, *xvss, *xdeff = NULL;
    int      *xier1, *xier2 = NULL;
    int      i, nnmax, mmj;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(m1);
    i3 = INTEGER_POINTER(m2);
    i4 = INTEGER_POINTER(m3);
    i5 = INTEGER_POINTER(m4);
    i6 = INTEGER_POINTER(iper);
    i7 = INTEGER_POINTER(jyear);
    i8 = INTEGER_POINTER(month);
    d2 = NUMERIC_POINTER(tau2);
    i9= INTEGER_POINTER(ns);
    i10 = INTEGER_POINTER(nfe);
    i11 = INTEGER_POINTER(npe);
    d3 = NUMERIC_POINTER(ar);
    i12 = INTEGER_POINTER(logt);
    i13 = INTEGER_POINTER(iopt);
    d4 = NUMERIC_POINTER(omin);
    d5 = NUMERIC_POINTER(omax);
    i14 = INTEGER_POINTER(nmax);
    i15 = INTEGER_POINTER(mj);

    nnmax = *i14;
    mmj = *i15;
    PROTECT(ans = allocVector(VECSXP, 8));
    SET_VECTOR_ELT(ans, 0, ff = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 1, over = allocVector(REALSXP, 1)); 
    SET_VECTOR_ELT(ans, 2, aic = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 3, xss = allocVector(REALSXP, nnmax*mmj)); 
    SET_VECTOR_ELT(ans, 4, vss = allocVector(REALSXP, nnmax*mmj*mmj)); 
    SET_VECTOR_ELT(ans, 5, deff = allocVector(REALSXP, nnmax)); 
    SET_VECTOR_ELT(ans, 6, ier1 = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 7, ier2 = allocVector(INTSXP, 1));

    d6 = NUMERIC_POINTER(ff);
    d7 = NUMERIC_POINTER(over);
    d8 = NUMERIC_POINTER(aic);
    d9 = NUMERIC_POINTER(xss);
    d10 = NUMERIC_POINTER(vss);
    d11 = NUMERIC_POINTER(deff);
    i16 = INTEGER_POINTER(ier1);
    i17 = INTEGER_POINTER(ier2);

    F77_CALL(seasonf) (d1,i1,i2,i3,i4,i5,i6,i7,i8,d2,i9,i10,i11,d3,i12,i13,d4,d5,i14,i15,d6,d7,d8,d9,d10,d11,i16,i17);

    xff = REAL(ff);
    xover = REAL(over);
    xaic = REAL(aic);
    xxss = REAL(xss);
    xvss = REAL(vss);
    xdeff = REAL(deff);
    xier1 = INTEGER(ier1);
    xier2 = INTEGER(ier2);

    *xff = *d6;
    *xover = *d7;
    *xaic = *d8;
    for(i=0; i<nnmax*mmj; i++) xxss[i] = d9[i];
    for(i=0; i<nnmax*mmj*mmj; i++) xvss[i] = d10[i];
    for(i=0; i<nnmax; i++) xdeff[i] = d11[i];
    *xier1 = *i16;
    *xier2 = *i17;

    UNPROTECT(1);

    return ans;
}

