#-----------------------------------------------------------------------------#
#                                                                             #
#            Interlaboratory Study Program ILS IN R                           #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores Sanchez                                       #
#              Professor of Mathematics Department                            #
#              Escuela Politecnica Nacional, Ecuador                          #
#              miguel.flores@epn.edu.ec                                       #
#                                                                             #
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Main function to create a 'ils.fqcdata' object
#-----------------------------------------------------------------------------#
##' Functional Quality Control Data
##' 
##' It Creates an object of class 'ils.fqcd' to perform statistical quality control. 
##' This object is used to plot functional data.
##' @aliases ils.fqcdata 
##' @param x A (n x m) matrix or data-frame. 
##' The m is the number of points observed in each curve, 
##' and n is the number of curves for each laboratory. 
##' @param p the number of laboratories.
##' @param index.laboratory is the laboratory index. The index.laboratory length should be equal a p 
##' @param argvals  Argvals, by default: 1:m.
##' @param rangeval Range of discretization points, by default: range(argvals).
##' @param names (optional) list with tree components: main an overall title, 
##' xlab title for x axis and ylab title for y axis.
##' @export
##' @references 
##' \describe{
##'   \item{}{Febrero-Bande, M. and Oviedo, M. (2012),
##'    "Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28.}
##'   \item{}{Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014), 
##'   "Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.}
##' }
##' @examples

##' library(ILS)
##' data(TG)
##' delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
##' fqcdata <- ils.fqcdata(TG, p = 7, argvals = delta)
##' xlab <- "Temperature (C)"
##' ylab <- "Mass (%)"
##' main <- "TG curves obtained from calcium oxalate"
##' plot(x = fqcdata, main = main, xlab=xlab, ylab=xlab,col = 1:7,legend = TRUE)

ils.fqcdata <- function(x, p = NULL, index.laboratory = NULL, argvals = NULL, rangeval = NULL, 
                        names = NULL)
#.........................................................................
{
  
  if (!is.data.frame(x) & !is.matrix(x) & is.null(x))
    stop("object must be a matrix or data.frame")

  if (is.null(p))
    stop("The number of laboratories should more than one")

  if (p == 1)
    stop("The number of laboratories should more than one")
  
  if (nrow(x) %% p != 0)
  stop("The number of laboratories must be multiple of number the rows")
  m <- dim(x)[2] # quality characteristics
  n <- dim(x)[1]/p # number of samples or observations
  
  if (is.null(argvals)) argvals <- seq(from = 1 ,to = m ,length.out = m )
  if (is.null(rangeval)) rangeval <- c(min(argvals),max(argvals))
  
  if (!is.null(index.laboratory )) 
    {if (length(index.laboratory) != p)
       stop("The index.laboratory length should be equal a p")}
  else 
  { index.laboratory <- paste("Lab",1:p)}
  
  ils.fdata <- fdata(mdata = x, argvals = argvals, rangeval = rangeval, names = names)
  
  
  oldClass(x) <- class(x) 

  result <- list (curves = as.matrix(x), ils.fdata = ils.fdata, 
                  index.laboratory = index.laboratory, p = as.numeric(p), n = n, m = m)
  
  oldClass(result) <- c("ils.fqcdata")
  
   
  return(result)
  
} # ils.fqcdata
#.........................................................................
