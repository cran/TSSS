# TSSS 1.3.4

* Removed C wrapper functions and registered entry points for the routines accessed by the `.Fortran` interface to call Fortran subroutines.

* URL: http://www.data.jma.go.jp/obd/stats/etrn/ moved to https://www.data.jma.go.jp/obd/stats/etrn/ in Rainfall.Rd and Temperature.Rd.

* Added a new argument `log.base = c("e" (default), "10")` to `season()`.

* Added error handling when local span is too small compared to max.arorder to lsar().

* Added a ‘`NEWS.md`’ file.


# TSSS 1.3.3

* Added a new argument `period` (period of one cycle) to `lsqr()`.

* Added a new argument `cov` (covariance matrix of the smoother) to `trend()`.

* Added a new argument `llkhood` (log-likelihood) to `ngsmth()`. 

* Fixed the number of frequencies for spectrum plots in `lsar()`.

* Fixed a bug in `tsmooth()` where the value of `esterr` was not included in the return value.

* Changed grayscale gradient in `pfilter()` and `pfilterNL()` plots.


# TSSS 1.3.2

* Added examples of `pairs` plots to the help pages for `HAKUSAN` and `Haibara` datasets.

* Fixed a bug in `season()` and changed the interpretation of argument `period`.

  In the case of `period = 7`, the data is daily (7 days a week) and not weekly.

* Corrected the source code of the Fortran subroutines FFARMA (armaftf.f) and FFSEAS (seasonf.f) to avoid runtime errors.

* Added error handling for missing data values to `arfit()`, `armafit()` and `armafit2()`.


# TSSS 1.3.1

* Exported functions (i.e., `print.tvvar()`, `plot.tvvar()` and `plot.pfilter()`) for an extension package **TSSSomp** using OpenMP.

* Removed the perspective plot from the output of `ngsmth()`.

* Added the following references:

  * Kitagawa, G. (2020)
   Introduction to Time Series Modeling with Applications in R.
   Chapman & Hall/CRC.

  * Kitagawa, G. (2020)
   Introduction to Time Series Modeling with R (in Japanese).
   Iwanami Publishing Company.


# TSSS 1.3.0

* Added a new argument `lag` to `period()`.

* Added `pfilter()` (Particle filtering and smoothing) and `pfilterNL()` (Particle filtering and nonlinear smoothing) functions.

* Added `armafit2()` function (Scalar ARMA model fitting: Estimate all ARMA models within the user-specified maximum order using maximum likelihood method).

* Added tsp attribute to all datasets.

* Added `Haibara`, `Nikkei225` and `Rainfall` datasets.

* Added an example using `Haibara` dataset in `crscor()`.

* Output the value of `arcoef` for each order of the AR model in `arfit()`.

* Changed value of `noisew` of examples in `ngsmth()`.

* Specified the allowed values for `trend.order`, `seasonal.order`, `ar.order` and `period` arguments in `season()`.

* Added `Note` about frequency setting, i.e., for data with the tsp attribute, the `frequency` attribute has priority over the `period` argument in `season()`.

* Changed the default value for `lag` from 10 to sqrt(n), where n is the length of the input time series.

* Changed the plot style of evolutionary power spectra in `plot.tvspc()`.

* Changed the function name from `armaimp()` to `armachar()`.

* Fixed a bug in subroutine istat3 (armafitf.f).

* Changed default values for `wminmax` in case of `noisew = 2` (Pearson) and `vminmax` in case of `noisev = 2` (Pearson) in `ngsim()`.


# TSSS 1.2.4

* Used the gammafn function provided by R instead of the `DGAMMA` function in the original Fortran source, because an error was found in TSSS-00install.out on r-patched-solaris-x86.


# TSSS 1.2.3

* Changed to use the `warning()` function instead of the `cat()` function to print warning messages.


# TSSS 1.2.2

* Fixed Fortran code for subroutine ARYULE (comsub.f) according to the warning message when using gfortran with `-fpic -g -O2 -Wall -mtune=native`.


# TSSS 1.2.1

* Changed the implementation of the random number generator to Mersenne-Twister.


# TSSS 1.2.0

* Fixed bugs in `arfit()`, `ngsmth()` and `trend()`.

* Fixed Fortran code according to the warning message when using gfortran with `-Wall -pedantic`.


* Merged `gdensty()`, `cdensty()`, `pdensty()`, `exdensty()`, `chi2densty()`, `dexdensty()` and `udensty()` into `pdfunc()`.

* Renamed `lsar1()` and `lsar2()` functions to `lsar()` and `lsar.chgpt()` respectively.

* Added S3 methods for class objects returned by functions.

* Added ‘`src/init.c`’.


# TSSS 1.1.0

* Changed `Sunspot` dataset into a part of `sunspot.year` in R`s core **datasets** package.

* Corrected the tile of the `HAKUSAN` dataset.


# TSSS 1.0.1

* Corrected ‘`inst/doc/index.html`’ according to the error message from the W3C Markup Validator Service.
