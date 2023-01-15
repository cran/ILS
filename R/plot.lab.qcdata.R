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

#-------------------------------------------------------------------------
# plot.lab.qcdata
#-------------------------------------------------------------------------
##' Plot method for 'lab.qcdata' objects
##'
##' Generic function for plotting objects of 'lab.qcdata' class. Results of  univariate ILS studies are graphically shown.
##'
##' @method plot lab.qcdata
##' @param x  An object of class \code{lab.qcdata} (Univariate Quality Control Data).
##' @param xlab Title for the x axis.
##' @param ylab Title for the y axis.
##' @param col Color of type material, when there only one.
##' @param ...  Other arguments to be passed to or from methods.
##' @export
##' @examples
##' library(ILS)
##' data(Glucose)
##' Glucose.qcdata <- lab.qcdata(Glucose)
##' str(Glucose.qcdata)
##' plot(Glucose.qcdata)

plot.lab.qcdata <- function(x, xlab = NULL, ylab = NULL, col = "blue", ...)
  #..............................................................................
{
  if(!is.null(x) & !inherits(x, "lab.qcdata") & !is.data.frame(x))
    stop("x must be an objects of class (or extending) 'lab.qcdata'")

  oldpar <- par(mar = c(4, 3, 1, 1) + 0.1)
  data.name <- attributes(x)$data.name

  if (is.null(xlab)) ylab <- "Laboratory"
  if (is.null(ylab)) xlab <- data.name


  m <- length(unique(x$material))

  material <- x$material
  laboratory <- x$laboratory
  var <- x$x


  if (m > 1){

    print(dotplot(laboratory ~ var,
                  groups = material,horizontal = T,
                  key = simpleKey(levels(as.factor(material)), space = "right"),
                  xlab = xlab,
                  aspect=0.5, ylab = ylab))

  }else{

    print(dotplot(laboratory ~ var,
                  horizontal = T,
                  xlab = xlab,
                  aspect=0.5, ylab=ylab,col = col))

  }

  par(oldpar)


  #.........................................................................
} # plot.lab.qcdata
#...........................................................................
