#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(ngsimf)  (int*,int*,int*,int*,int*,int*,int*,int*,double*,double*,double*,int*,double*,double*,double*,int*,double*,double*,double*);

SEXP NgsimC(SEXP m1, SEXP m2, SEXP m3, SEXP m, SEXP k, SEXP n, SEXP ini, SEXP noisew, SEXP wmin, SEXP wmax, SEXP pw, SEXP noisev, SEXP vmin, SEXP vmax, SEXP pv, SEXP period, SEXP arcoef, SEXP x)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10;

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
    i8 = INTEGER_POINTER(noisew);
    d1 = NUMERIC_POINTER(wmin);
    d2 = NUMERIC_POINTER(wmax);
    d3 = NUMERIC_POINTER(pw);
    i9 = INTEGER_POINTER(noisev);
    d4 = NUMERIC_POINTER(vmin);
    d5 = NUMERIC_POINTER(vmax);
    d6 = NUMERIC_POINTER(pv);
    i10 = INTEGER_POINTER(period);
    d7 = NUMERIC_POINTER(arcoef);
    d8 = NUMERIC_POINTER(x);

    nn = *i6;
    PROTECT(ans = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(ans, 0, y = allocVector(REALSXP, nn));

    d9 = NUMERIC_POINTER(y);

    F77_CALL(ngsimf) (i1,i2,i3,i4,i5,i6,i7,i8,d1,d2,d3,i9,d4,d5,d6,i10,d7,d8,d9);

    yy = REAL(y);

    for(i=0; i<nn; i++) yy[i] = d9[i];
    UNPROTECT(1);

    return ans;
}

