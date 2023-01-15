#-----------------------------------------------------------------------------#
#                                                                             #
#            Interlaboratory Study Program ILS IN R                           #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores Sanchez                                      #
#              Professor of Mathematics Department                            #
#              Escuela Politecnica Nacional, Ecuador                          #
#              miguel.flores@epn.edu.ec                                       #
#                                                                             #
#-----------------------------------------------------------------------------#
#-------------------------------------------------------------------------
# plot.ils.fqcd
#-------------------------------------------------------------------------
##' Plotting method for 'ils.fqcdata' objects
##'
##' Generic function to plot objects of 'ils.fqcdata' class
##'
##' @method plot ils.fqcdata
##' @param x  An object class of\code{ils.fqcdata} (Functional Quality Control Data)
##' @param type 1-character string giving the type of plot desired.
##' The following values are possible for fdata class object: "l" for lines (by default),
##' "p" for points, , "o" for overplotted points and lines, "b", "c" for (empty if "c")
##' points joined by lines, "s" and "S" for stair steps and "h" for histogram-like
##' vertical lines. Finally, "n" does not produce any points or lines.
##' The following values are possible for fdata2d class object: "image.contour" (by default) to display three-dimensional data and add the contour lines, "image" to display three-dimensional data, "contour" to display a contour plot, "persp" to display a perspective plots of a surface over the x-y plane and "filled.contour" to display a contour plot with the areas between the contours filled in solid color.
##' @param main Main title for the plot.
##' @param xlab Title for the x axis.
##' @param ylab Title for the y axis.
##' @param ylim The y limits of the plot.
##' @param x.co It speficies the x co-ordinates to be used to place a legend.
##' @param y.co It specifies the y co-ordinates to be used to place a legend.
##' @param legend Logical argument. Default is TRUE then The legend default is used.
##' @param col Color specifications
##' @param ...  Other arguments passed to matplot function (for fdata class).
##' @export
##' @references
##' Febrero-Bande, M. and Oviedo, M. (2012),
##' "Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28.
##' Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014),
##' "Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.
##' 
##' @examples
##' \dontrun{
##' library(ILS)
##' data(TG)
##' delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
##' fqcdata <- ils.fqcdata(TG, p = 7, argvals = delta)
##' windows()
##' xlab <- "Temperature (C)"
##' ylab <- "Mass (%)"
##' main <- "TG curves obtained from calcium oxalate"
##' plot(x = fqcdata, main = main, xlab=xlab, ylab=ylab,legend = TRUE)}

plot.ils.fqcdata <- function(x,type = "l", main = NULL, xlab = NULL,
                          ylab = NULL, ylim = NULL,
                          x.co = NULL, y.co = NULL,
                          legend = TRUE, col = NULL, ...)
  #..............................................................................
{

  if(!is.null(x) & !inherits(x, "fdata") & !inherits(x, "ils.fqcdata"))
    stop("x must be an objects of class (or extending) ils.fqcdata")

  objectfdata <- x$ils.fdata

  data <- objectfdata[["data"]]
  tt <- argvals(objectfdata)
  ylim <- range(data) + c(-2,2)

  if (is.null(x.co)) x.co <- min(tt)
  if (is.null(y.co)) y.co <- 0.99 * max(data)

  p <- x$p
  n <- x$n

  if (is.null(col)) col <- terrain.colors(p)

  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)
  objectfdata$names <- list(main = main, xlab= xlab, ylab = ylab)
  plot.fdata(objectfdata,col = col, ylim = ylim, ...)

  rect(par("usr")[1],
       par("usr")[3],
       par("usr")[2],
       par("usr")[4],
       col  =  "white")
  box(col  =  "#CCCCCC")
  grid(col  =  "#CCCCCC")

  for(i in 0:(p-1))
    lines(objectfdata[(i*n)+(1:n),],col=col[i+1],...)

  if (legend == TRUE)
  legend(x.co, y.co, unique(x$index.laboratory), col = col, lty = 1, lwd = 2, cex =0.8)
  par(oldpar)
  #.........................................................................
} # plot.ilsfqcdata
#...........................................................................
