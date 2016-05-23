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
# plot.lab.fqcd
#-------------------------------------------------------------------------
##' Plotting method for 'lab.fqcd' objects
##' 
##' Generic function to plot objects of 'lab.fqcd' class
##' 
##' @method plot lab.fqcd
##' @param x  Object lab.fqcd (Functional Quality Control Data)
##' @param main Main title for the plot
##' @param xlab Title for the x axis
##' @param ylab Title for the y axis
##' @param ylim The y limits of the plot
##' @param x.co It speficies the x co-ordinates to be used to place a legend.
##' @param y.co It specifies the y co-ordinates to be used to place a legend.
##' @param legend Logical argument. Default is TRUE then The legend default is used. 
##' @param col Color specifications
##' @param ...  Arguments to be passed to or from methods.
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
##' curves.fqcd <- lab.fqcd(TG, argvals = delta)
##' xlab <- "Temperature/ C"
##' ylab <- "Mass/ %"
##' main <- "TG curves obtained from calcium oxalate"
##' p <- dim(curves.fqcd$curves)[3]
##' col <- 1:p
##' plot(x = curves.fqcd, main, xlab, ylab, col= col,legend = FALSE)
##' legend(45,70,paste("Lab",col),
##'       col = col, lty = 1, lwd = 2, cex =0.8)
##'

plot.lab.fqcd <- function(x, main = NULL, xlab = NULL, ylab = NULL, ylim = NULL, x.co = NULL, y.co = NULL, 
                          legend = TRUE, col = NULL, ...)
  #..............................................................................      
{

  if(!is.null(x) & !inherits(x, "fdata") & class(x)!="lab.fqcd")
    stop("x must be an objects of class (or extending) lab.fqcd")
  
  if (is.null(ylim)) ylim <- range(x$curves) + c(-2,2)
  if (is.null(main)) main <- attributes(x)$type.data
  if (is.null(xlab)) xlab <- "argvals"
  if (is.null(ylab)) ylab <- "Curves"
  
  data <- x$curves.fdata[[1]]$data
  tt <- x$curves.fdata[[1]]$argvals
  
  if (is.null(x.co)) x.co <- min(tt)
  if (is.null(y.co)) y.co <- 0.99 * max(data)
  
  p <- length(x$curves.fdata)
  if (is.null(col)) col <- terrain.colors(p)
  
  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)
  plot.fdata(x[[2]][[1]], xlab = xlab, ylab = ylab, col = 1, type = "l",
             cex.lab = 0.7, cex.axis = 0.7, main = main, ylim = ylim,...)
  
  rect(par("usr")[1],
       par("usr")[3],
       par("usr")[2],
       par("usr")[4],
       col  =  "white")
  box(col  =  "#CCCCCC")
  grid(col  =  "#CCCCCC")
  
  for(i in 1:p)
    lines(x[[2]][[i]],col=col[i],...)
  
  if (legend==TRUE)
  legend(x.co, y.co, paste("Lab",1:p),col = col, lty = 1, lwd = 2, cex =0.8) 
  par(oldpar)
  #.........................................................................
} # plot.list.fdata
#...........................................................................

#-------------------------------------------------------------------------
# plot.lab.fqcs
#-------------------------------------------------------------------------
##' Plotting method for 'lab.fqcs' objects
##' 
##' Generic function to plot objects of 'lab.fqcs' class. Results of functional ILS studies are graphically shown.
##' 
##' @method plot lab.fqcs
##' @param x  Object functional data or a list with objects of functional data type 
##' @param main Main title for the plot
##' @param xlab Title for the x axis
##' @param ylab Title for the y axis
##' @param ylim The y limits of the plot
##' @param x.co It speficies the x co-ordinates to be used to place a legend.
##' @param y.co It specifies the y co-ordinates to be used to place a legend.
##' @param legend Logical argument. Default is TRUE then The legend default is used. 
##' @param col Color specifications
##' @param ...  arguments to be passed to or from methods
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
##' curves.fqcd <- lab.fqcd(TG, argvals = delta)
##' curves.fqcs <- lab.fqcs(curves.fqcd)
##' summary(curves.fqcs)
##' names(curves.fqcs)
##' class(curves.fqcs$mean.i)
##' xlab <- "Temperature/ C"
##' ylab <- "Mass/ %"
##' main <- "Functional Mean Estimation by Laboratory"
##' p <- dim(curves.fqcd$curves)[3]
##' col <- 1:p
##' plot(curves.fqcs$mean.i,main = main, xlab = xlab, ylab = ylab, col = col,legend = FALSE)
##' legend(45,70,paste("Lab",1:p),
##'       col = col,lty = 1,lwd = 2,cex = 0.7)
plot.lab.fqcs <- function(x, main = NULL, xlab = NULL, ylab = NULL, ylim = NULL, x.co = NULL, y.co = NULL,
                          legend = TRUE, col = NULL, ...)
  #..............................................................................      
{
    if(!is.null(x) & !inherits(x, "fdata") & !inherits(x, "h.fqcs") & 
       !inherits(x, "k.fqcs") & !inherits(x, "lab.fqcs" ) & !inherits(x, "w.fqcs" ) & !inherits(x, "d.fqcs" ))
    stop("x must be an objects of class (or extending) fqcs")

  type.data <- attributes(x)$object.name
  
  if (is.null(main)) main <- "Statistical Functional"
  if (is.null(xlab)) xlab <- "argvals"
  if (is.null(ylab)) ylab <- "Curves"

  if (type.data!="fdata"){
    p <- x$p
    if (is.null(col)) col <- terrain.colors(p)
    }
  
  result <- switch(type.data, 
                  "data" =  {plot.lab.fqcd(x,main = main,xlab = xlab,ylab = ylab, ylim = ylim,
                                           x.co = x.co, y.co=y.co, legend = legend, col = col,...)} ,
                  "h.fqcs" =   {plot.lab.fqcd(x$hi,main = main,xlab = xlab,ylab = ylab, ylim = ylim,
                                              x.co = x.co, y.co=y.co, legend = legend, col = col,...)},
                  "k.fqcs" ={plot.lab.fqcd(x$ki,main = main,xlab = xlab,ylab = ylab, ylim = ylim,
                                           x.co = x.co, y.co=y.co, legend = legend, col = col,...)},
                  "QuantileDepth" ={
                    quantile <- lab.fqcd(rbind(x["low"]$low$data,x["up"]$up$data))
                    plot.lab.fqcd(quantile,main = main,xlab = xlab,ylab = ylab, ylim = ylim,
                                  x.co = x.co, y.co=y.co, legend = legend, col = col,...)},
                  "QuantileWalter" ={
                    if (length(x)==2){
                      quantile <- lab.fqcd(rbind(x["low"]$low$data,x["up"]$up$data))
                    }else{
                      quantile <- lab.fqcd(x["up"]$up$data)
                    }
                    plot.lab.fqcd(quantile,main = main,xlab = xlab,ylab = ylab, ylim = ylim,
                                  x.co = x.co, y.co=y.co, legend = legend, col = col,...)},
                  "lab.fqcs" ={ 
                    
                    par(mfrow=c(2,2))
                    
                    plot.lab.fqcd(x$mean.i,main = main,xlab = xlab,ylab = ylab,ylim = ylim,
                                  x.co = x.co, y.co=y.co, legend = legend, col = col, 
                                  sub = "Functional Mean Estimation by Laboratory",...)
                    
                    plot.lab.fqcd(x$s2.i,main = main,xlab = xlab,ylab = ylab,ylim = ylim,
                                  x.co = x.co, y.co=y.co, col = col, legend = FALSE,
                                  sub = "Functional Variance Estimation by Laboratory",...)
  
                    plot.fdata(x$mean, xlab = xlab, ylab = ylab, col = "red", type = "l",
                               cex.lab = 0.7, cex.axis = 0.7, main = main, pch  =  16, 
                               axes  =  TRUE,sub = "Functional Mean Global Estimation",...)
                    rect(par("usr")[1],
                         par("usr")[3],
                         par("usr")[2],
                         par("usr")[4],
                         col  =  "white")
                    box(col  =  "#CCCCCC")
                    grid(col  =  "#CCCCCC")
                    lines(x$mean, col = col,...)

                    plot.fdata(x$S2, xlab = xlab, ylab = ylab, col = "red", type = "l",
                               cex.lab = 0.7, cex.axis = 0.7, main = main, pch  =  16, 
                               axes  =  TRUE, sub = "Functional Variance  Global Estimation",...)
                    rect(par("usr")[1],
                         par("usr")[3],
                         par("usr")[2],
                         par("usr")[4],
                         col  =  "white")
                    box(col  =  "#CCCCCC")
                    grid(col  =  "#CCCCCC")
                    lines(x$S2, col = col,...)
                    par(mfrow=c(1,1))
                  })
  #.........................................................................
} # plot.lab.fqcs
#...........................................................................

