#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(marspcf)(int*, int*, double*, double*, int*, Rcomplex*, double*, double*, double*, double*, double*);

SEXP MarspcC(SEXP m, SEXP l, SEXP a, SEXP e, SEXP nf)
{
    double *d1, *d2, *d3, *d4, *d5, *d6, *d7;
    int      *i1, *i2, *i3;
    Rcomplex *c1;

    SEXP ans = R_NilValue, p = R_NilValue, amp = R_NilValue, ang = R_NilValue, coh = R_NilValue, fnc = R_NilValue, frnc = R_NilValue;
    Rcomplex *xp = NULL;
    double *xamp, *xang, *xcoh, *xfnc, *xfrnc = NULL;
    int  nf1, ll, llnf1, i;

    i1 = INTEGER_POINTER(m);
    i2 = INTEGER_POINTER(l);
    d1 = NUMERIC_POINTER(a);
    d2 = NUMERIC_POINTER(e);
    i3 = INTEGER_POINTER(nf);

    ll = (*i2) * (*i2);
    nf1 = *i3 + 1;
    llnf1 = ll * nf1;
    PROTECT(ans = allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, p = allocVector(CPLXSXP, llnf1));
    SET_VECTOR_ELT(ans, 1, amp = allocVector(REALSXP, llnf1));
    SET_VECTOR_ELT(ans, 2, ang = allocVector(REALSXP, llnf1));
    SET_VECTOR_ELT(ans, 3, coh = allocVector(REALSXP, llnf1));
    SET_VECTOR_ELT(ans, 4, fnc = allocVector(REALSXP, llnf1));
    SET_VECTOR_ELT(ans, 5, frnc = allocVector(REALSXP, llnf1));

    xp = COMPLEX(p);
    xamp = REAL(amp);
    xang = REAL(ang);
    xcoh = REAL(coh);
    xfnc = REAL(fnc);
    xfrnc = REAL(frnc);

    c1 = COMPLEX_POINTER(p);
    d3 = NUMERIC_POINTER(amp);
    d4 = NUMERIC_POINTER(ang);
    d5 = NUMERIC_POINTER(coh);
    d6 = NUMERIC_POINTER(fnc);
    d7 = NUMERIC_POINTER(frnc);

    F77_CALL(marspcf) (i1,i2,d1,d2,i3,c1,d3,d4,d5,d6,d7);

    for(i=0; i<llnf1; i++) xp[i] = c1[i];
    for(i=0; i<llnf1; i++) xamp[i] = d3[i];
    for(i=0; i<llnf1; i++) xang[i] = d4[i];
    for(i=0; i<llnf1; i++) xcoh[i] = d5[i];
    for(i=0; i<llnf1; i++) xfnc[i] = d6[i];
    for(i=0; i<llnf1; i++) xfrnc[i] = d7[i];

    UNPROTECT(1);

    return ans;
}
