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
# plot.ils.fqcs
#-------------------------------------------------------------------------
##' Plotting method for 'ils.fqcs' objects
##'
##' Generic function to plot objects of 'ils.fqcs' class. Results of functional ILS studies are graphically shown.
##'
##' @method plot ils.fqcs
##' @param x  An object of class \code{ils.fqcs} (Functional Quality Control Statistics).
##' @param type 1-character string giving the type of plot desired.
##' The following values are possible for fdata class object: "l" for lines (by default),
##' "p" for points, , "o" for overplotted points and lines, "b", "c" for (empty if "c")
##' points joined by lines, "s" and "S" for stair steps and "h" for histogram-like
##' vertical lines. Finally, "n" does not produce any points or lines.
##' The following values are possible for fdata2d class object: "image.contour" (by default) to display three-dimensional data and add the contour lines, "image" to display three-dimensional data, "contour" to display a contour plot, "persp" to display a perspective plots of a surface over the x-y plane and "filled.contour" to display a contour plot with the areas between the contours filled in solid color.
##' @param xlab Title for the x axis.
##' @param ylab Title for the y axis.
##' @param legend Logical argument. Default is TRUE then The legend default is used.
##' @param col Color specifications.
##' @param ...  Other arguments to be passed to or from methods.
##' @export
##' @references
##' Febrero-Bande, M. and Oviedo, M. (2012),
##' "Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28.
##' 
##' Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014),
##'   "Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.
##' 
##' @examples
##' library(ILS)
##' data(TG)
##' delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
##' fqcdata <- ils.fqcdata(TG, p = 7, argvals = delta)
##' xlab <- "Temperature/ C"
##' ylab <- "Mass/ %"
##' fqcstat <- ils.fqcs(fqcdata)
##' plot(fqcstat, xlab = xlab, ylab = ylab,legend = TRUE)


plot.ils.fqcs <- function(x, type = "l", xlab = NULL, ylab = NULL,
                          legend = TRUE, col = NULL,  ...)
  #..............................................................................
{
    if(!is.null(x) & !inherits(x, "ils.fqcs" ))
    stop("x must be an objects of class (or extending) ils.fqcs")

  p <- x$p

  if (is.null(col)) col <- terrain.colors(p)


                    par(mfrow=c(2,2))

                    data <- x$mean.i$data
                    tt <- x$mean.i$argvals
                    main <- "Functional Mean by Laboratory"
                    x.co <- min(tt)
                    y.co <- 0.99 * max(data)
                    x$mean.i$names <- list(main = main, xlab= xlab, ylab = ylab)

                    oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)
                    plot.fdata(x$mean.i, col = col, type = type,
                               cex.lab = 0.7, cex.axis = 0.7,...)

                    rect(par("usr")[1],
                         par("usr")[3],
                         par("usr")[2],
                         par("usr")[4],
                         col  =  "white")
                    box(col  =  "#CCCCCC")
                    grid(col  =  "#CCCCCC")

                    lines(x$mean.i,col = col,...)

                  if (legend == TRUE)
                    legend(x.co, y.co, unique(x$fqcdata$index.laboratory),col = col,
                           lty = 1, lwd = 2, cex =0.8)
                  par(oldpar)

                  ############ S2.i

                  data <- x$S2.i$data
                  tt <- x$S2.i$argvals
                  main <- "Functional Variance by Laboratory"
                  x.co <- min(tt)
                  y.co <- 0.99 * max(data)
                  x$S2.i$names <- list(main = main, xlab= xlab, ylab = ylab)

                  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)
                  plot.fdata(x$S2.i, xlab = xlab, ylab = ylab, col = col, type = type,
                             cex.lab = 0.7, cex.axis = 0.7, main = main,...)

                  rect(par("usr")[1],
                       par("usr")[3],
                       par("usr")[2],
                       par("usr")[4],
                       col  =  "white")
                  box(col  =  "#CCCCCC")
                  grid(col  =  "#CCCCCC")

                  lines(x$S2.i,col = col,...)

                  if (legend == TRUE)
                    legend(x.co, y.co, unique(x$fqcdata$index.laboratory),col = col, lty = 1, lwd = 2, cex =0.8)
                  par(oldpar)

                  ###### mean

                  data <- x$mean$data
                  tt <- x$mean$argvals
                  main <- "Global Functional Mean"
                  x.co <- min(tt)
                  y.co <- 0.99 * max(data)
                  x$mean$names <- list(main = main, xlab= xlab, ylab = ylab)

                  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)
                  plot.fdata(x$mean, xlab = xlab, ylab = ylab, col = col, type = "l",
                             cex.lab = 0.7, cex.axis = 0.7, main = main, ...)

                  rect(par("usr")[1],
                       par("usr")[3],
                       par("usr")[2],
                       par("usr")[4],
                       col  =  "white")
                  box(col  =  "#CCCCCC")
                  grid(col  =  "#CCCCCC")

                  lines(x$mean,col = col,...)

                  par(oldpar)
                  ######### S2

                  data <- x$S2$data
                  tt <- x$S2$argvals
                  main <- "Global Functional Variance"
                  x.co <- min(tt)
                  y.co <- 0.99 * max(data)
                  x$S2$names <- list(main = main, xlab= xlab, ylab = ylab)

                  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)
                  plot.fdata(x$S2, xlab = xlab, ylab = ylab, col = col, type = "l",
                             cex.lab = 0.7, cex.axis = 0.7, main = main, ...)

                  rect(par("usr")[1],
                       par("usr")[3],
                       par("usr")[2],
                       par("usr")[4],
                       col  =  "white")
                  box(col  =  "#CCCCCC")
                  grid(col  =  "#CCCCCC")

                  lines(x$S2,col = col,...)

                  par(oldpar)

                  par(mfrow=c(1,1))

  #.........................................................................
} # plot.ils.fqcs
#...........................................................................


#-------------------------------------------------------------------------
# plot.mandel.fqcs
#-------------------------------------------------------------------------
##' Plotting method for 'mandel.fqcs' objects
##'
##' Generic function to plot objects of 'mandel.fqcs' class. Results of functional ILS studies are graphically shown.
##'
##' @method plot mandel.fqcs
##' @param x  An object of class \code{mandel.fqcs}.
##' @param xlab Title for the x axis.
##' @param ylab Title for the y axis.
##' @param x.co It speficies the x co-ordinates to be used to place a legend.
##' @param y.co It specifies the y co-ordinates to be used to place a legend.
##' @param legend Logical argument. Default is TRUE then The legend default is used.
##' @param col Color specifications.
##' @param ...  Other arguments to be passed to or from methods.
##' @export
##' @references
##' Febrero-Bande, M. and Oviedo, M. (2012),
##' "Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28.
##' 
##' Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014),
##' "Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.
##' 
##' @examples
##' \dontrun{
##' library(ILS)
##' data(TG)
##' delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
##' fqcdata <- ils.fqcdata(TG, p = 7, argvals = delta)
##' mandel.tg <- mandel.fqcs(fqcdata.tg,nb = 200)
##' plot(mandel.tg,legend = F,col=c(rep(3,5),1,1))
##'   }

plot.mandel.fqcs <- function(x, xlab = NULL, ylab = NULL,x.co = NULL, y.co =NULL,
                          legend = TRUE, col = NULL, ...)
  #..............................................................................
{
  if(!is.null(x) & !inherits(x, "mandel.fqcs"))
    stop("x must be an objects of class (or extending) mandel.fqcs")
  est <- x
  fdata.m <- est$m.ils.fdata
  fqcs <- est$fqcs
  index <- fqcs$fqcdata$index.laboratory
  p <- fqcs$fqcdata$p
  n <- fqcs$fqcdata$n
  mdata <- fqcs$fqcdata$m

  #statistics
  S2 <- fqcs$S2
  mean <- fqcs$mean
  mean.i <- fqcs$mean.i
  S2.i <- fqcs$S2.i
  Sbar <- fqcs$Sbar
  S2r <- fqcs$S2r
  ils.h <- fqcs$ils.h
  ils.k <- fqcs$ils.k
  n.h <- fqcs$d_h.m
  n.k <- fqcs$d_k.m

  distboot.k <- est$distboot.k
  dist.k <- est$dist.k
  fmax.k <- est$fmax.k
  fmin.k <- est$fmin.k

  distboot.h <- est$distboot.h
  dist.h <- est$dist.h
  fmax.h <- est$fmax.h
  fmin.h <- est$fmin.h


  #' # $H(t)$ y $K(t)$
  #'
  ## ------------------------------------------------------------------------


  par(mfrow=c(2,2))
  ## dH
  main <- expression(d_H)
  if (is.null(xlab)) xlab <- "Laboratory"
  if (is.null(ylab)) ylab <- "Statistic"
  if (is.null(col)) col <- terrain.colors(p)

  ylim <- c(0,max(n.h,dist.h) + 2)
  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)


  barplot(height = t(n.h),width = 0.05,main=main,names.arg = index,
          beside = T,
          col=col,xlab=xlab,
          ylim=ylim,axisnames=T)

  abline(h = dist.h, col = "red")

  par(oldpar)

  ### H
  data <- est$m.fest.h$data
  tt <- est$m.fest.h$argvals
  main <- "H(x) statistic"
  if (is.null(x.co)) x.co <- min(tt)
  if (is.null(y.co)) y.co <- 0.99 * 3
  ylim <-  c(-3,3)

  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)

  plot(est$m.fest.h,xlab = xlab, ylab = ylab, col = col, type = "l", ylim = ylim,
       cex.lab = 0.7, cex.axis = 0.7, main = main,...)

  rect(par("usr")[1],
       par("usr")[3],
       par("usr")[2],
       par("usr")[4],
       col  =  "white")
  box(col  =  "#CCCCCC")
  grid(col  =  "#CCCCCC")

  lines(est$m.fest.h,col =col)
  lines(fmax.h, lwd=1,col="red")
  lines(fmin.h, lwd=1,col="red")

  if (legend == TRUE)
    legend(x.co, y.co, c(index,"CL H(x)"),
                         col = c(col,"red"), lty = 1, lwd = 2, cex =0.8)
  par(oldpar)

  ### dK
  main <- expression(d_K)

  ylim <- c(0,max(n.k,dist.k) + 2)
  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)


  barplot(height = t(n.k),width = 0.05,main=main,names.arg = index,
          beside = T,
          col=col,xlab=xlab,
          ylim=ylim,axisnames=T)

  abline(h = dist.k, col = "red")

  par(oldpar)

  ### K
  data <- est$m.fest.k$data
  tt <- est$m.fest.k$argvals
  main <- "K(x) statistic"
  x.co <- min(tt)
  y.co <- 0.99 * 3

  ylim <-  c(0,3)

  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)

  plot(est$m.fest.k,xlab = xlab, ylab = ylab, col = col, type = "l", ylim = ylim,
       cex.lab = 0.7, cex.axis = 0.7, main = main,...)

  rect(par("usr")[1],
       par("usr")[3],
       par("usr")[2],
       par("usr")[4],
       col  =  "white")
  box(col  =  "#CCCCCC")
  grid(col  =  "#CCCCCC")

  lines(est$m.fest.k,col =col)
  lines(fmax.k, lwd=1,col="red")

  if (legend == TRUE)
    legend(x.co, y.co, c(index,"CL K(x)"),
           col = c(col,"red"), lty = 1, lwd = 2, cex =0.8)
  par(oldpar)


  par(mfrow=c(1,1))
  #.........................................................................
} # plot.mandel.fqcs
#...........................................................................


