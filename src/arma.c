#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(armaf)(int*, int*, double*, double*, double*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, int*, int*);

SEXP arma(SEXP arorder, SEXP maorder, SEXP arcoef, SEXP macoef, SEXP v, SEXP n, SEXP lag, SEXP kmax, SEXP nf)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;

    SEXP ans = R_NilValue, g = R_NilValue, acov = R_NilValue, parcor = R_NilValue, spec = R_NilValue;
    SEXP roota = R_NilValue, rootb = R_NilValue, ier = R_NilValue, jer = R_NilValue;
    double *xg, *xacov, *xparcor, *xspec, *xroota, *xrootb = NULL;
    int  *nier, *njer = NULL;
    int ar, ma, lag0, lag1, kmax1, nf1, ar2, ma2;
    int i;

    i1 = INTEGER_POINTER(arorder);
    i2 = INTEGER_POINTER(maorder);
    d1 = NUMERIC_POINTER(arcoef);
    d2 = NUMERIC_POINTER(macoef);
    d3 = NUMERIC_POINTER(v);
    i3 = INTEGER_POINTER(n);
    i4 = INTEGER_POINTER(lag);
    i5 = INTEGER_POINTER(kmax);
    i6 = INTEGER_POINTER(nf);

    ar = *i1;
    ma = *i2;
    lag0 = *i4;
    lag1 = lag0 + 1;
    kmax1 = *i5+1;
    nf1 = *i6 + 1;
    ar2 = ar*2;
    if( ar2 == 0 ) ar2 = 2;
    ma2 = ma*2;
    if( ma2 == 0 ) ma2 = 2;

    PROTECT(ans = allocVector(VECSXP, 8));
    SET_VECTOR_ELT(ans, 0, g = allocVector(REALSXP, kmax1));
    SET_VECTOR_ELT(ans, 1, acov = allocVector(REALSXP, lag1));
    SET_VECTOR_ELT(ans, 2, parcor = allocVector(REALSXP, lag0));
    SET_VECTOR_ELT(ans, 3, spec = allocVector(REALSXP, nf1));
    SET_VECTOR_ELT(ans, 4, roota = allocVector(REALSXP, ar2));
    SET_VECTOR_ELT(ans, 5, rootb = allocVector(REALSXP, ma2));
    SET_VECTOR_ELT(ans, 6, ier = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 7, jer = allocVector(INTSXP, 1));

    d4 = NUMERIC_POINTER(g);
    d5 = NUMERIC_POINTER(acov);
    d6 = NUMERIC_POINTER(parcor);
    d7 = NUMERIC_POINTER(spec);
    d8 = NUMERIC_POINTER(roota);
    d9 = NUMERIC_POINTER(rootb);
    i7 = INTEGER_POINTER(ier);
    i8 = INTEGER_POINTER(jer);

    F77_CALL(armaf) (i1,i2,d1,d2,d3,i3,i4,i5,i6,d4,d5,d6,d7,d8,d9,i7,i8);

    xg = REAL(g);
    xacov = REAL(acov);
    xparcor = REAL(parcor);
    xspec = REAL(spec);
    xroota = REAL(roota);
    xrootb = REAL(rootb);
    nier = INTEGER(ier);
    njer = INTEGER(jer);

    for(i=0; i<kmax1; i++) xg[i] = d4[i];
    for(i=0; i<lag1; i++) xacov[i] = d5[i];
    for(i=0; i<lag0; i++) xparcor[i] = d6[i];
    for(i=0; i<nf1; i++) xspec[i] = d7[i];
    if( ar == 0 ) xroota = NULL;
    if( ar != 0 ) for(i=0; i<ar2; i++)  xroota[i] = d8[i];
    if( ma == 0 ) xrootb = NULL;
    if( ma != 0 ) for(i=0; i<ma2; i++)  xrootb[i] = d9[i];
    *nier = *i7;
    *njer = *i8;

    UNPROTECT(1);

    return ans;
}
