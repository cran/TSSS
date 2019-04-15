#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(armasp) (double*,int*,double*,int*,double*,int*,double*);

SEXP arsp(SEXP a, SEXP m, SEXP b, SEXP l, SEXP sig2, SEXP nf)
{
    double *d1,*d2,*d3,*d4;
    int *i1,*i2,*i3;

    SEXP ans =  R_NilValue, spec = R_NilValue;

    double *xspec = NULL;
    int    i, nf1;

    d1 = NUMERIC_POINTER(a);
    i1 = INTEGER_POINTER(m);
    d2 = NUMERIC_POINTER(b);
    i2 = INTEGER_POINTER(l);
    d3 = NUMERIC_POINTER(sig2);
    i3 = INTEGER_POINTER(nf);

    nf1 = *i3 + 1;
    PROTECT(ans = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(ans, 0, spec = allocVector(REALSXP, nf1)); 

    d4 = NUMERIC_POINTER(spec);

    F77_CALL(armasp) (d1,i1,d2,i2,d3,i3,d4);

    xspec = REAL(spec);

    for(i=0; i<nf1; i++) xspec[i] = d4[i];

    UNPROTECT(1);

    return ans;
}
