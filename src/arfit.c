#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(arfitf)(double*, int*, int*, int*, int*, int*, double*, double*, int*, double*, double*, double* );

SEXP ArfitC(SEXP y, SEXP n, SEXP lag, SEXP nf, SEXP mj2, SEXP isw)
{
    double *d1,*d2,*d3,*d4,*d5,*d6;
    int *i1,*i2,*i3,*i4,*i5,*i6;

    SEXP ans = R_NilValue, sig = R_NilValue, aic = R_NilValue, mar = R_NilValue, arcoef = R_NilValue, parcor = R_NilValue, spec = R_NilValue;
    double *xsig, *xaic, *xarcoef, *xparcor, *xspec = NULL;
    int  *xmar = NULL;
    int i, lg, lg1, nf1;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(lag);
    i3 = INTEGER_POINTER(nf);
    i4 = INTEGER_POINTER(mj2);
	i5 = INTEGER_POINTER(isw);

    lg = *i2;
    lg1 = lg + 1;
    nf1 = *i3 + 1;

    PROTECT(ans = allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, sig = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 1, aic = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 2, mar = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 3, arcoef = allocVector(REALSXP, lg*lg));
    SET_VECTOR_ELT(ans, 4, parcor = allocVector(REALSXP, lg));
    SET_VECTOR_ELT(ans, 5, spec = allocVector(REALSXP, nf1));

    d2 = NUMERIC_POINTER(sig);
    d3 = NUMERIC_POINTER(aic);
    i6 = INTEGER_POINTER(mar);
    d4 = NUMERIC_POINTER(arcoef);
    d5 = NUMERIC_POINTER(parcor);
    d6 = NUMERIC_POINTER(spec);

    F77_CALL(arfitf) (d1,i1,i2,i3,i4,i5,d2,d3,i6,d4,d5,d6);

    xsig = REAL(sig);
    xaic = REAL(aic);
    xmar = INTEGER(mar);
    xarcoef = REAL(arcoef);
    xparcor = REAL(parcor);
    xspec = REAL(spec);

    for(i=0; i<lg1; i++) xsig[i] = d2[i];
    for(i=0; i<lg1; i++) xaic[i] = d3[i];
    *xmar = *i6;
    for(i=0; i<lg*lg; i++) xarcoef[i] = d4[i];
    for(i=0; i<lg; i++) xparcor[i] = d5[i];
    for(i=0; i<nf1; i++) xspec[i] = d6[i];

    UNPROTECT(1);

    return ans;
}
