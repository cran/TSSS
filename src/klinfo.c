#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(klinfof) ( int*, double*, int*, double*, double*, double*, int*, double*, double*, double* );

SEXP KlinfoC(SEXP distg, SEXP paramg, SEXP distf, SEXP paramf, SEXP xmin, SEXP xmax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7;
    int   *i1,*i2,*i3;

    SEXP ans = R_NilValue, nint = R_NilValue, dx = R_NilValue, fkli = R_NilValue, gint = R_NilValue;
    double *xdx, *xfkli, *xgint = NULL;
    int *xnint = NULL;
    int   i;

    i1 = INTEGER_POINTER(distg);
    d1 = NUMERIC_POINTER(paramg);
    i2 = INTEGER_POINTER(distf);
    d2 = NUMERIC_POINTER(paramf);
    d3 = NUMERIC_POINTER(xmin);
    d4 = NUMERIC_POINTER(xmax);

    PROTECT(ans = allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ans, 0, nint = allocVector(INTSXP, 4));
    SET_VECTOR_ELT(ans, 1, dx = allocVector(REALSXP, 4));
    SET_VECTOR_ELT(ans, 2, fkli = allocVector(REALSXP, 4));
    SET_VECTOR_ELT(ans, 3, gint = allocVector(REALSXP, 4));

    i3 = INTEGER_POINTER(nint);
    d5 = NUMERIC_POINTER(dx);
    d6 = NUMERIC_POINTER(fkli);
    d7 = NUMERIC_POINTER(gint);

    F77_CALL(klinfof) (i1,d1,i2,d2,d3,d4,i3,d5,d6,d7);

    xnint = INTEGER(nint);
    xdx = REAL(dx);
    xfkli = REAL(fkli);
    xgint = REAL(gint);

    for(i=0; i<4; i++) xnint[i] = i3[i];
    for(i=0; i<4; i++) xdx[i] = d5[i];
    for(i=0; i<4; i++) xfkli[i] = d6[i];
    for(i=0; i<4; i++) xgint[i] = d7[i];

    UNPROTECT(1);

    return ans;
}
