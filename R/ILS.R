##' Interlaboratory Study
##'
##' It performs interlaboratory studies (ILS) to detect those laboratories that provide non-consistent results when comparing to others.
##' It permits to work simultaneously with various testing materials, from standard univariate, and functional data analysis (FDA) perspectives.
##' The univariate approach based on ASTM E691-08 consist of estimating the Mandel's h and k statistics to identify those laboratories
##' that provide more significant different results, testing also the presence of outliers by Cochran and Grubbs tests, Analysis of variance (ANOVA)
##' techniques are provided (F and Tuckey tests) to test differences in the testing variable means corresponding to test differences in means corresponding to differente laboratories per each material.
##' Taking into account the functional nature of data retrieved in analytical chemistry, applied physics and engineering (spectra, thermograms, etc.).
##' ILS package provides a FDA approach for  functional Mandel's k and h statistics by  smoothing bootstrap resampling of distribution.
##'
##' @name ILS
##'
##' @aliases ILS
##'
##' @docType package
##'
##' @title Interlaboratoty Study
##'
##' @import multcomp
##' @import depthTools
##' @import fda.usc
##' @import MASS
##' @import graphics
##' @importFrom grDevices terrain.colors
##' @importFrom utils setTxtProgressBar str txtProgressBar
##' @importFrom stats aov confint pf pt qf qt sd var
NULL

##' @title Glucose in Serum
##'
##' @description
##' Dataset corresponding to serum glucose (measurements of glucose concentration in blood used to control diabetes) testing.
##' Eight laboratories conducted tests to five different blood samples tagged with different references, ranging them from low  sugar content to very high.
##' Three replicates were obtained for each sample. It is retrieved from ASTM E 691 standard.
##'
##' @name Glucose
##'
##' @docType data
##
##' @format A data frame with 120 observations composed of the following 4 variables:
##' \describe{
##'   \item{Glucose}{Glucose content in Serum}
##'   \item{Replicate}{Number of glucose measurement corresponding to each material}
##'   \item{Material}{Level of glucose, ranging from low content of sugar to very high level of glucose in blood.}
##'   \item{Laboratory}{Laboratories conducted tests}
##' }
##' @keywords datasets
##'
##' @references
##' \describe{
##'   \item{}{ASTM E 691 (1999). Standard  practice  for  conducting  an  interlaboratory  study  to  determine
##' the precision of a test method. American Society for Testing and Materials. West Conshohocken, PA, USA. }
##'}
##' @examples
##' library(ILS)
##' data(Glucose)
##' summary(Glucose)
##' attach(Glucose)
##' str(Glucose)
##' table(Replicate,Material,Laboratory)
##' table(Laboratory,Material)
##' st <- with(Glucose, tapply(Glucose, list(Material,Laboratory), mean))
##' st
##'
NULL

##' @title Dataset composed of the initial decomposition temperature (IDT) of different samples of Calcium Oxalate, obtained by 7 different laboratories
##'
##' @description Initial decomposition temperature (IDT) is a parameter defined by temperature at which a material loss 5\% of its weight
##' when it is heated using a constant rate. One hundred and five calcium oxalate samples were tested by thermogravimetric analysis (TG),
##' obtaining 105 TG curves from which the IDT is extracted. Summarizing, IDT dataset is composed of the IDT values of
##' calcium oxalate obtained by 7 different laboratories that analyze 15 oxalate samples each one.:
##' Laboratory 1 uses a simultaneous thermal analyzer (STA) with an old calibration program, Laboratory 2 to Laboratory 4 use a SDT simultaneous analyzer,
##' Laboratory 6 utilizes a SDT simultaneous analyzer with an old calibration, and Laboratory 7 uses a SDT simultaneous analyzer with a biased calibration
##' (2 degrees Celsius shifted from the zinc melting point).
##'
##' @name IDT
##'
##' @docType data
##'
##' @format Dataframe of dimension 105 x 44. The first column corresponds to IDT variable, the second (Sample) is the replicate number,
##'  the third is the tested material (Material), and fourth is the laboratory.
##'
##' @references
##' \describe{
##'   \item{}{Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014),
##'   "Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.}
##' }
##' @examples
##' library(ILS)
##' data(IDT)
##' summary(IDT)
##' attach(IDT)
##' str(IDT)
##' table(Sample,Run,Laboratory)
##' table(Laboratory,Run)
##' st <- with(IDT, tapply(IDT, list(Run,Laboratory), mean))
##' st

NULL


##' @title Thermogravimetry curves
##'
##' @description One hundred and five Calcium oxalate samples were tested by thermogravimetric (TG) analysis,
##' obtaining 105 TG curves that shows the mass loss of oxalate depending on time when samples are heated at a constant temperature rate.
##' Dataset is composed by fifteen TG curves of 1000 observations each of overall 7 different laboratories.
##' Laboratory 1 uses a simultaneous thermal analyzer (STA) with an old calibration program, Laboratory 2 to Laboratory 4
##' use a SDT simultaneous analyzer, Laboratory 6 utilizes a SDT simultaneous analyzer with an old calibration, and Laboratory 7
##'  uses a SDT simultaneous analyzer with a biased calibration (2 degrees Celsius shifted from the zinc melting point).
##'
##' @name TG
##'
##' @docType data
##'
##' @format A 15 x 1000 x 7 dimension array ,  where each matrix consists of the 15 TG curves obtained testing 15 different oxalate samples,
##' and evaluated in 1000 different values of temperature. These 15 curves were obtained for each of the overall 7 laboratories that have performed the esperiments.
##'
##' @references
##' \describe{
##'   \item{}{Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014),
##'   "Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.}
##' }
##' @examples
##' library(ILS)
##' data(TG)
##' summary(TG)
##'
NULL
