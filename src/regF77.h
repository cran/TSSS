/*
*  TSSS : R Package for Time Series Analysis with State Space Model
*  Copyright (C) 2013    The Institute of Statistical Mathematics
*
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program; if not, write to the Free Software
*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
*
*
*  ismrp at grp.ism.ac.jp
*/

#include <R.h>
#include <Rinternals.h>
#include <libintl.h>

#define _(String) (String)

/* Fortran : */

void F77_NAME(arfit)(double *y, int *n, int *lag, int *nf, int *mj2, int *isw,
              double *sig2, double *aic, int *mar, double *a, double *par,
			  double *sp);

void F77_NAME(arma)(int *arorder, int *maorder, double *arcoef, double *macoef,
              double *v, int *lag, int *kmax, int *nf, double *g, double *acov,
			  double *parcor, double *spec, double *roota, double *rootb,
			  int *ier, int *jer);
	
void F77_NAME(armaft)(double *y, int *n, int *m, int *l, int *mlmax,
              int *iparam, double *ar0, double *ma0, double *sig2, double *flk,
			  double *aic, double *ar, double *ma, int *ier);

void F77_NAME(armaft2)(double *y, int *n, int *mmax, int *lmax, int *mlmax,
              double *sig2, double *flk, double *aic, double *ar, double *ma,
			  int *ier);

void F77_NAME(boxcoxf)(double *y, int *n, double *aiczt, double *ffzt,
              double *aicz, double *ffz, double *mean, double *var, double *zz);

void F77_NAME(crscorf)(double *y, int *n, int *l, int *lag, double *outmin,
              double *outmax, double *cov, double *cor, double *mean);

void F77_NAME(klinfof)(int *distg, double *paramg, int *distf, double *paramf,
              double *xmin, double *xmax, int *nint, double *dx, double *fkli,
			  double *gint);

void F77_NAME(lsar1)(double *y, int *n, int *lag, int *ns0, int *nb, int *nf,
              int *ns, int *n0, int *n1, int *iflag, int *ms, double *sds,
			  double *aics, int *mp, double *sdp, double *aicp, double *as,
			  int *mfs, double *sig2s, int *nnf, int *ier);

void F77_NAME(lsar2)(double *y, int *n, int *k, int *n0, int *n1, int *n2,
              int *ne, double *aics, double *aicmin, int *nmin);

void F77_NAME(lsqr)(double *y, int *n, int *k, int *period, int *mj,
              double *aic, double *sig2, int *imin, double *reg, double *data);

void F77_NAME(marfit)(double *y, int *n, int *l, int *lag,
              double *amin, double *vmin, double *aic, int *mmin);

void F77_NAME(marlsq)(double *y, int *n, int *l, int *lag,
              double *arcoef, double *v, int *lmax, double *aic);

void F77_NAME(marspcf)(int *m, int *l, double *a, double *e, int *nf,
              Rcomplex *p, double *amp, double *ang, double *coh, double *fnc,
			  double *frnc);

void F77_NAME(ngsmthf)(double *y, int *n, int *noisev, double *tau2, double *bv,
              int *noisew, double *sig2, double *bw, int *initd, double *trend,
			  double *smt, double *lkhood, int *ns, int *nfe, int *npe, int *k1);

void F77_NAME(denstyf)(int *model, double *param, double *xmin, double *xmax,
              int *k, double *f);

void F77_NAME(periodf)(double *y, int *n, int *np, int *iw, int *lag,
              double *outmin, double *outmax, double *pe, double *spe, int *ifg);

void F77_NAME(fftperf)(double *y, int *n, int *iw,
              double *pe, double *spe, int *np, int *ifg);

void F77_NAME(pfilter)(double *y, int *n, int *m, int *model, int *lag,
              int *inid, double *sig2, double *tau2, double *alpha,
			  double *bigtau2, double *sig2i, double *xmin, double *xmax,
			  int *ix, double *t, double *ff);

void F77_NAME(pfiltern)(double *y, int *n, int *m, int *lag, double *sig2,
              double *tau2, double *xmin, double *xmax, int *ix, double *t,
			  double *ff);

void F77_NAME(polreg)(double *y, int *n, int *k,
              double *a, double *sig2, double *aic, int *imin, double *data);

void F77_NAME(season)(double *y, int *n, int *m1, int *m2, int *m3, int *m4,
              int *iper, int *jyear, int *month, double *tau2, int *ns,
			  int *nfe, int *npe, double *ar, int *logt, int *iopt,
			  double *omin, double *omax, int *nmax, int *mj, double *ff,
			  double *over, double *aic, double *xss, double *vss, double *deff,
			  int *ier1, int *ier2);

void F77_NAME(simssmf)(int *m1, int *m2, int *m3, int *m, int *k, int *n,
              int *ini, double *sig2, int *period, double *tau1, double *tau2,
			  double *tau3, double *arcoef, double *x, double *y);

void F77_NAME(ngsimf)(int *m1, int *m2, int *m3, int *m, int *k, int *n,
              int *ini, int *noisew, double *wmin, double *wmax, double *pw,
			  int *noisev, double *vmin, double *vmax, double *pv, int *period,
			  double *arcoef, double *x, double *y);

void F77_NAME(trend)(double *y, int *n, int *m, int *iopt, double *tau20,
              double *delta, double *tau2, double *sig2, double *lkhood,
			  double *aic, double *xss, double *vss, double *rs);

void F77_NAME(smoothf)(double *y, int *n, int *m, int *k, double *f, double *g,
              double *h, double *q, double *r, double *xf, double *vf, int *nfe,
			  int *npe, double *omin, double *omax, int *nmiss, int *n0,
			  int *nn, double *xss, double *vss, double *lkhood, double *aic);
	
void F77_NAME(tvvarf)(double *y, int *n, int *m, double *tau20, int *iopt,
              double *delta, double *tvvar, double *nordat, double *y1, int *n1,
			  double *trend, double *noise, double *taumax, double *sig2m,
			  double *ffmax, double *aic);

void F77_NAME(tvar)(double *y, int *n, int *arorder, int *torder, int *span,
              int *method, int *nout, int *outlier, double *tau20,
			  double *delta, double *tau2, double *sigma2, double *lkhood,
			  double *aic, double *arcoef, double *parcor);

void F77_NAME(tvspc)(int *n, int *arorder, int *span, int *nf, int *ivar,
              double *sigma2, double *arcoef, double *var, double *spec);

void F77_NAME(unicorf)(double *y, int *n, int *lag, double *outmin,
              double *outmax, double *cov, double *mean);

void F77_NAME(armasp)(double *a, int *m, double *, int *l, double *sig2,
              int *nf, double *spec);
