#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(pfilterf)(double*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, int*, double*, double*);

SEXP PfiltC(SEXP y, SEXP n, SEXP m, SEXP model, SEXP lag, SEXP inid, SEXP sig2, SEXP tau2, SEXP alpha, SEXP bigtau2, SEXP sig2i, SEXP xmin, SEXP xmax, SEXP ix)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10;
    int *i1,*i2,*i3,*i4,*i5,*i6;

    SEXP ans = R_NilValue, t = R_NilValue, ff = R_NilValue;
    double *xt, *xff = NULL;
    int   i, nn;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(m);
    i3 = INTEGER_POINTER(model);
    i4 = INTEGER_POINTER(lag);
    i5 = INTEGER_POINTER(inid);
    d2 = NUMERIC_POINTER(sig2);
    d3 = NUMERIC_POINTER(tau2);
    d4 = NUMERIC_POINTER(alpha);
    d5 = NUMERIC_POINTER(bigtau2);
    d6 = NUMERIC_POINTER(sig2i);
    d7 = NUMERIC_POINTER(xmin);
    d8 = NUMERIC_POINTER(xmax);
    i6 = INTEGER_POINTER(ix);

    nn = *i1;
    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, t = allocVector(REALSXP, 8*nn));
    SET_VECTOR_ELT(ans, 1, ff = allocVector(REALSXP, 1));

    d9 = NUMERIC_POINTER(t);
    d10 = NUMERIC_POINTER(ff);

    F77_CALL(pfilterf) (d1,i1,i2,i3,i4,i5,d2,d3,d4,d5,d6,d7,d8,i6,d9,d10);

    xt = REAL(t);
    xff = REAL(ff);

    for(i=0; i<8*nn; i++) xt[i] = d9[i];
    *xff = *d10;

    UNPROTECT(1);

    return ans;
}

