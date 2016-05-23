#-----------------------------------------------------------------------------#
#                                                                             #
#            Interlaboratory Study Program ILS IN R                           #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores Sanchez                                      #
#              Professor of Mathematics Department                            #
#              Escuela Politecnica Nacional, Ecuador                           #
#              miguel.flores@epn.edu.ec                                       #
#                                                                             #
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Main function to create a 'lab.fqcd' object
#-----------------------------------------------------------------------------#
##' Functional Quality Control Data
##' 
##' It Creates an object of class 'lab.fqcd' to perform statistical quality control. 
##' This object is used to plot functional data.
##' @aliases lab.fqcd 
##' @param x A (m x p) matrix or data-frame. Alternatively an (m x p x n) array. 
##' The m parameter is the number of curves, p defines the number of points observed in each curve, 
##' and n is the number of replicates. 
##' @param argvals  Argvals, by default: 1:p.
##' @param rangeval Range of discretization points, by default: range(argvals).
##' @export
##' @references 
##' \describe{
##'   \item{}{Febrero-Bande, M. and Oviedo, M. (2012),
##'    "Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28.}
##'   \item{}{Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014), 
##'   "Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.}
##' }
##' @examples
##' \dontrun{
##' library(ILS)
##' data(TG) 
##' delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
##' curves.fqcd <- lab.fqcd(TG, argvals = delta)
##' windows()
##' xlab <- "Temperature (0C)"
##' ylab <- "Mass (%)"
##' main <- "TG curves obtained from calcium oxalate"
##' p <- curves.fqcd$p
##' col <- terrain.colors(p)
##' plot(x = curves.fqcd, main, xlab, ylab,legend = FALSE,col = col)
##' legend(45,70,c(paste("Lab",c(1:7))),
##'       col = col,lty = 1, lwd = c(rep(1,7),2), cex = 0.7)
##'       }
lab.fqcd <- function(x, argvals = NULL, rangeval = NULL)
#.........................................................................
{
  if (!is.array(x) & !is.data.frame(x) & !is.matrix(x))
    stop("object must be a array or a matrix or data.frame")

  
  curves.fdata <- list()
  
  if (length(dim(x))==3){
    m <- dim(x)[2] # quality characteristics
    n <- dim(x)[1] # number of samples or observations
    p <- dim(x)[3]       # # number of laboratories
    oldClass(x) <- c("array")
  } else{
    m <- dim(x)[2] # quality characteristics
    n <- 1 # number of samples or observations
    p <- dim(x)[1]       # number of laboratories
    oldClass(x) <- c("matrix")
  }

  if (is.null(argvals)) argvals <- seq(from = 1 ,to = m ,length.out = m )
  if (is.null(rangeval)) rangeval <- c(min(argvals),max(argvals))
  
  for (i in 1:p){
    if (length(dim(x))==3){
       curves.fdata[[i]] <- fdata(mdata = x[,,i], argvals = argvals, rangeval = rangeval)
    }else{
      curves.fdata[[i]] <- fdata(mdata = x[i,], argvals = argvals, rangeval = rangeval)
    }
}
  
  oldClass(curves.fdata) <- c("list.fdata") 
  result <- list (curves = x, curves.fdata = curves.fdata, p = p, n = n, m = m)
  
  attr(result, "object.name") <-"data" 
  attr(result, "type.data") <- "Functional Data List"
  
  oldClass(result) <- c("lab.fqcd")
  
   
  return(result)
  
} # lab.fdata
#.........................................................................
