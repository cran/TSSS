#include <R.h>
#include <Rdefines.h>
#include "TSSS.h"

extern void F77_NAME(lsar1f) (double*,int*,int*,int*,int*,int*,
int*,int*,int*,int*,int*,double*,double*,int*,double*,double*,double*,int*,double*,int*);

SEXP Lsar1C(SEXP y, SEXP n, SEXP lag, SEXP ns0, SEXP nb, SEXP nf)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13;

    SEXP ans =  R_NilValue, ns = R_NilValue, n0 = R_NilValue, n1 = R_NilValue, iflag = R_NilValue, ms = R_NilValue, sds = R_NilValue, aics = R_NilValue, mp = R_NilValue, sdp = R_NilValue, aicp = R_NilValue, as = R_NilValue, mfs = R_NilValue, sig2s = R_NilValue, nnf = R_NilValue;

    double *xsds, *xaics, *xsdp, *xaicp, *xas, *xsig2s = NULL;
    int    *xns, *xn0, *xn1, *xiflag, *xms, *xmp, *xmfs, *xnnf = NULL;
    int    i, lg, nnb;

    d1 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(n);
    i2 = INTEGER_POINTER(lag);
    i3 = INTEGER_POINTER(ns0);
    i4 = INTEGER_POINTER(nb);
    i5 = INTEGER_POINTER(nf);

    lg = *i2;
    nnb = *i4;
    PROTECT(ans = allocVector(VECSXP, 14));
    SET_VECTOR_ELT(ans, 0, ns = allocVector(INTSXP, nnb));
    SET_VECTOR_ELT(ans, 1, n0 = allocVector(INTSXP, nnb));
    SET_VECTOR_ELT(ans, 2, n1 = allocVector(INTSXP, nnb));
    SET_VECTOR_ELT(ans, 3, iflag = allocVector(INTSXP, nnb));
    SET_VECTOR_ELT(ans, 4, ms = allocVector(INTSXP, nnb));
    SET_VECTOR_ELT(ans, 5, sds = allocVector(REALSXP, nnb)); 
    SET_VECTOR_ELT(ans, 6, aics = allocVector(REALSXP, nnb));
    SET_VECTOR_ELT(ans, 7, mp = allocVector(INTSXP, nnb));
    SET_VECTOR_ELT(ans, 8, sdp = allocVector(REALSXP, nnb)); 
    SET_VECTOR_ELT(ans, 9, aicp = allocVector(REALSXP, nnb));
    SET_VECTOR_ELT(ans, 10, as = allocVector(REALSXP, nnb*lg));
    SET_VECTOR_ELT(ans, 11, mfs = allocVector(INTSXP, nnb));
    SET_VECTOR_ELT(ans, 12, sig2s = allocVector(REALSXP, nnb));
    SET_VECTOR_ELT(ans, 13, nnf = allocVector(INTSXP, nnb));

    i6 = INTEGER_POINTER(ns);
    i7 = INTEGER_POINTER(n0);
    i8 = INTEGER_POINTER(n1);
    i9 = INTEGER_POINTER(iflag);
    i10 = INTEGER_POINTER(ms);;
    d2 = NUMERIC_POINTER(sds);
    d3 = NUMERIC_POINTER(aics);
    i11 = INTEGER_POINTER(mp);
    d4 = NUMERIC_POINTER(sdp);
    d5 = NUMERIC_POINTER(aicp);
    d6 = NUMERIC_POINTER(as);
    i12 = INTEGER_POINTER(mfs);
    d7 = NUMERIC_POINTER(sig2s);
    i13 = INTEGER_POINTER(nnf);

    F77_CALL(lsar1f) (d1,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,d2,d3,i11,d4,d5,d6,i12,d7,i13);

    xns = INTEGER(ns);
    xn0 = INTEGER(n0);
    xn1 = INTEGER(n1);
    xiflag = INTEGER(iflag);
    xms = INTEGER(ms);
    xsds = REAL(sds);
    xaics = REAL(aics);
    xmp = INTEGER(mp);
    xsdp = REAL(sdp);
    xaicp = REAL(aicp);
    xas = REAL(as);
    xmfs = INTEGER(mfs);
    xsig2s = REAL(sig2s);
    xnnf = INTEGER(nnf);

    for(i=0; i<nnb; i++) xns[i] = i6[i];
    for(i=0; i<nnb; i++) xn0[i] = i7[i];
    for(i=0; i<nnb; i++) xn1[i] = i8[i];
    for(i=0; i<nnb; i++) xiflag[i] = i9[i];
    for(i=0; i<nnb; i++) xms[i] = i10[i];
    for(i=0; i<nnb; i++) xsds[i] = d2[i];
    for(i=0; i<nnb; i++) xaics[i] = d3[i];
    for(i=0; i<nnb; i++) xmp[i] = i11[i];
    for(i=0; i<nnb; i++) xsdp[i] = d4[i];
    for(i=0; i<nnb; i++) xaicp[i] = d5[i];
    for(i=0; i<lg*nnb; i++) xas[i] = d6[i];
    for(i=0; i<nnb; i++) xmfs[i] = i12[i];
    for(i=0; i<nnb; i++) xsig2s[i] = d7[i];
    for(i=0; i<nnb; i++) xnnf[i] = i13[i];

    UNPROTECT(1);

    return ans;
}
