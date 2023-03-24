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

#include "regF77.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* .Fortran calls */

static const R_FortranMethodDef FortEntries[] = {
    {"arfit",    (DL_FUNC) &F77_NAME(arfit),    12},
    {"arma",     (DL_FUNC) &F77_NAME(arma),     16},
    {"armaft",   (DL_FUNC) &F77_NAME(armaft),   14},
    {"armaft2",  (DL_FUNC) &F77_NAME(armaft2),  11},
    {"boxcoxf",  (DL_FUNC) &F77_NAME(boxcoxf),   9},
    {"crscorf",  (DL_FUNC) &F77_NAME(crscorf),   9},
    {"klinfof",  (DL_FUNC) &F77_NAME(klinfof),  10},
    {"lsar1",    (DL_FUNC) &F77_NAME(lsar1),    21},
    {"lsar2",    (DL_FUNC) &F77_NAME(lsar2),    10},
    {"lsqr",     (DL_FUNC) &F77_NAME(lsqr),     10},
    {"marfit",   (DL_FUNC) &F77_NAME(marfit),    8},
    {"marlsq",   (DL_FUNC) &F77_NAME(marlsq),    8},
    {"marspcf",  (DL_FUNC) &F77_NAME(marspcf),  11},
    {"ngsmthf",  (DL_FUNC) &F77_NAME(ngsmthf),  16},
    {"denstyf",  (DL_FUNC) &F77_NAME(denstyf),   6},
    {"periodf",  (DL_FUNC) &F77_NAME(periodf),  10},
    {"fftperf",  (DL_FUNC) &F77_NAME(fftperf),   7},
    {"pfilter",  (DL_FUNC) &F77_NAME(pfilter),  16},
    {"pfiltern", (DL_FUNC) &F77_NAME(pfiltern), 11},
    {"polreg",   (DL_FUNC) &F77_NAME(polreg),    8},
    {"season",   (DL_FUNC) &F77_NAME(season),   28},
    {"simssmf",  (DL_FUNC) &F77_NAME(simssmf),  15},
    {"ngsimf",   (DL_FUNC) &F77_NAME(ngsimf),   19},
    {"trend",    (DL_FUNC) &F77_NAME(trend),    13},
    {"smoothf",  (DL_FUNC) &F77_NAME(smoothf),  22},
    {"tvvarf",   (DL_FUNC) &F77_NAME(tvvarf),   16},
    {"tvar",     (DL_FUNC) &F77_NAME(tvar),     16},
    {"tvspc",    (DL_FUNC) &F77_NAME(tvspc),     9},
    {"unicorf",  (DL_FUNC) &F77_NAME(unicorf),   7},
    {"armasp",   (DL_FUNC) &F77_NAME(armasp),    7},
    {NULL, NULL, 0}
};

void attribute_visible R_init_TSSS(DllInfo *dll) 
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
