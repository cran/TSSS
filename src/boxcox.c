#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(boxcoxf) ( double*, int*, double*, double*, double*, double*, double*, double*, double* );

SEXP BoxcoxC(SEXP y, SEXP n)
{
    double *d1,*d2,*d3,*d4,*d5, *d6, *d7, *d8;
    int   *i1;

    SEXP ans =  R_NilValue, aiczt = R_NilValue, ffzt = R_NilValue, aicz = R_NilValue, ffz = R_NilValue, mean = R_NilValue, var = R_NilValue, z = R_NilValue;
    double *xaiczt, *xffzt, *xaicz, *xffz, *xmean, *xvar, *xz = NULL;
    int   i, nn;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);

    nn = *i1;
    PROTECT(ans = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(ans, 0, aiczt = allocVector(REALSXP, 21));
    SET_VECTOR_ELT(ans, 1, ffzt = allocVector(REALSXP, 21));
    SET_VECTOR_ELT(ans, 2, aicz = allocVector(REALSXP, 21));
    SET_VECTOR_ELT(ans, 3, ffz = allocVector(REALSXP, 21));
    SET_VECTOR_ELT(ans, 4, mean = allocVector(REALSXP, 21));
    SET_VECTOR_ELT(ans, 5, var = allocVector(REALSXP, 21));
    SET_VECTOR_ELT(ans, 6, z = allocVector(REALSXP, nn));

    d2 = NUMERIC_POINTER(aiczt);
    d3 = NUMERIC_POINTER(ffzt);
    d4 = NUMERIC_POINTER(aicz);
    d5 = NUMERIC_POINTER(ffz);
    d6 = NUMERIC_POINTER(mean);
    d7 = NUMERIC_POINTER(var);
    d8 = NUMERIC_POINTER(z);

    F77_CALL(boxcoxf) (d1,i1,d2,d3,d4,d5,d6,d7,d8);

    xaiczt = REAL(aiczt);
    xffzt = REAL(ffzt);
    xaicz = REAL(aicz);
    xffz = REAL(ffz);
    xmean = REAL(mean);
    xvar = REAL(var);
    xz = REAL(z);

    for(i=0; i<21; i++) xaiczt[i] = d2[i];
    for(i=0; i<21; i++) xffzt[i] = d3[i];
    for(i=0; i<21; i++) xaicz[i] = d4[i];
    for(i=0; i<21; i++) xffz[i] = d5[i];
    for(i=0; i<21; i++) xmean[i] = d6[i];
    for(i=0; i<21; i++) xvar[i] = d7[i];
    for(i=0; i<nn; i++) xz[i] = d8[i];

    UNPROTECT(1);

    return ans;
}
