% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ils.fqcs.R
\name{ils.fqcs}
\alias{ils.fqcs}
\alias{ils.fqcs.default}
\alias{ils.fqcs.ils.fqcdata}
\alias{print.ils.fqcs}
\alias{summary.ils.fqcs}
\title{It developes an object of class 'ils.fqcs'}
\usage{
ils.fqcs(x, ...)

\method{ils.fqcs}{default}(x, argvals = NULL, rangeval = NULL, ...)

\method{ils.fqcs}{ils.fqcdata}(x, ...)

\method{print}{ils.fqcs}(x, ...)

\method{summary}{ils.fqcs}(object, ...)
}
\arguments{
\item{x}{An object of class \code{ils.fqcs} for which a print is desired.}

\item{...}{Other arguments passed to or from methods.}

\item{argvals}{Argvals, by default: \code{1:m}.}

\item{rangeval}{The range of discretization points, by default: range(argvals).}

\item{object}{An object of class \code{ils.fqcs} for which a summary is desired.}
}
\description{
Create an object of class 'ils.fqcs' to perform statistical quality control.
This function is used to compute requested FDA.
}
\examples{
library(ILS)
data(TG)
delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
fqcdata <- ils.fqcdata(TG, p = 7, argvals = delta)
xlab <- "Temperature/ C"
ylab <- "Mass/ \%"
fqcstat <- ils.fqcs(fqcdata)
plot(fqcstat, xlab = xlab, ylab = ylab,legend = TRUE)
}
\references{
Febrero-Bande, M. and Oviedo, M. (2012),
"Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28

Cuevas A., Febrero-Bande, M. and Fraiman, R. (2006), "On the use of the bootstrap for estimating functions with functional data".
Computational Statistics & Data Analysis 51, 2, 1063-1074.

Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014),
"Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.
}
