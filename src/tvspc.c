#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(tvspcf) (int*, int*, int*, int*, int*, double*, double*, double*, double*);

SEXP TvspcC(SEXP n, SEXP aorder, SEXP span, SEXP nf, SEXP ivar, SEXP sigma2, SEXP arcoef, SEXP var)
{
    double *d1,*d2,*d3,*d4;
    int    *i1,*i2,*i3,*i4,*i5;

    SEXP ans =  R_NilValue,  spec = R_NilValue;
    double *xspec = NULL;
    int   i, k;

    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(aorder);
    i3 = INTEGER_POINTER(span);
    i4 = INTEGER_POINTER(nf);
    i5 = INTEGER_POINTER(ivar);
    d1 = NUMERIC_POINTER(sigma2);
    d2 = NUMERIC_POINTER(arcoef);
    d3 = NUMERIC_POINTER(var);

    k = (*i4 + 1) * (*i1);
    PROTECT(ans = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(ans, 0, spec = allocVector(REALSXP, k));

    d4 = NUMERIC_POINTER(spec);

    F77_CALL(tvspcf) (i1,i2,i3,i4,i5,d1,d2,d3,d4);

    xspec = REAL(spec);

    for(i=0; i<k; i++) xspec[i] = d4[i];

    UNPROTECT(1);

    return ans;
}
