% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ILS.R
\docType{package}
\name{ILS}
\alias{ILS}
\title{Interlaboratoty Study}
\description{
Interlaboratory Study
}
\details{
It performs interlaboratory studies (ILS) to detect those laboratories that provide non-consistent results when comparing to others.
It permits to work simultaneously with various testing materials, from standard univariate, and functional data analysis (FDA) perspectives.
The univariate approach based on ASTM E691-08 consist of estimating the Mandel's h and k statistics to identify those laboratories
that provide more significant different results, testing also the presence of outliers by Cochran and Grubbs tests, Analysis of variance (ANOVA)
techniques are provided (F and Tuckey tests) to test differences in the testing variable means corresponding to test differences in means corresponding to differente laboratories per each material.
Taking into account the functional nature of data retrieved in analytical chemistry, applied physics and engineering (spectra, thermograms, etc.).
ILS package provides a FDA approach for  functional Mandel's k and h statistics by  smoothing bootstrap resampling of distribution.
}
