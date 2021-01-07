#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(marfitf)(double*, int*, int*, int*, double*, double*, double*, int*);

SEXP MarfitC(SEXP y, SEXP n, SEXP l, SEXP lag)
{
    double *d1, *d2, *d3, *d4;
    int      *i1, *i2, *i3, *i4;

    SEXP ans = R_NilValue, amin = R_NilValue, vmin = R_NilValue, aic = R_NilValue, mmin = R_NilValue;
    double *xamin, *xvmin, *xaic = NULL;
    int  *xmmin = NULL;
    int  ll, lg, lg1, i;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(l);
    i3 = INTEGER_POINTER(lag);

    ll = (*i2) * (*i2);
    lg = *i3;
    lg1 = lg + 1;
    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, amin = allocVector(REALSXP, ll*lg));
    SET_VECTOR_ELT(ans, 1, vmin = allocVector(REALSXP, ll));
    SET_VECTOR_ELT(ans, 2, aic = allocVector(REALSXP, lg1));
    SET_VECTOR_ELT(ans, 3, mmin = allocVector(INTSXP, 1));

    d2 = NUMERIC_POINTER(amin);
    d3 = NUMERIC_POINTER(vmin);
    d4 = NUMERIC_POINTER(aic);
    i4 = INTEGER_POINTER(mmin);

    F77_CALL(marfitf) (d1,i1,i2,i3,d2,d3,d4,i4);

    xamin = REAL(amin);
    xvmin = REAL(vmin);
    xaic = REAL(aic);
    xmmin = INTEGER(mmin);

    for(i=0; i<ll*lg; i++) xamin[i] = d2[i];
    for(i=0; i<ll; i++) xvmin[i] = d3[i];
    for(i=0; i<lg1; i++) xaic[i] = d4[i];
    *xmmin = *i4;

    UNPROTECT(1);

    return ans;
}
