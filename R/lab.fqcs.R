#-----------------------------------------------------------------------------#
#                                                                             #
#            Interlaboratory Study Program ILS IN R                           #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores Sanchez                                       #
#              Professor of Mathematics Department                            #
#              Escuela Politecnica Nacional, Ecuador                           #
#              miguel.flores@epn.edu.ec                                       #
#                                                                             #
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Main function to create a 'lab.fqcs' object
#-----------------------------------------------------------------------------#
##' It developes an object of class 'lab.fqcs'
##' 
##' Create an object of class 'lab.fqcs' to perform statistical quality control. 
##' This function is used to compute requested FDA.
##'
## @param x  Object lab.fqcd (Functional Quality Control Data)
##' @export
##' @references
##' \describe{
##'   \item{}{Febrero-Bande, M. and Oviedo, M. (2012),
##'    "Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28.}
##'   \item{}{Cuevas A., Febrero-Bande, M. and Fraiman, R. (2006), "On the use of the bootstrap for estimating functions with functional data". 
##'   Computational Statistics & Data Analysis 51, 2, 1063-1074. }
##'   \item{}{Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014), 
##'   "Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.}
##' }
##' @examples
##' \dontrun{
##' library(ILS)
##' data(TG)
##' delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
##' curves.fqcd <- lab.fqcd(TG, argvals = delta)
##' p <- curves.fqcd$p
##' curves.fqcs <- lab.fqcs(curves.fqcd)
##' summary(curves.fqcs)
##' names(curves.fqcs)
##' 
##' ##### Statistical Functional
##' windows(20,10)
##' xlab <- "Temperature (C)"
##' ylab <- "Mass (%)"
##' main <- "Statistical Functional"
##' col <- 1:p
##' }
lab.fqcs <- function(x, ...) {
  UseMethod("lab.fqcs")
}

##' @rdname lab.fqcs
##' @method lab.fqcs default
##' @inheritParams lab.fqcd
##' @param ... Arguments passed to or from methods.
lab.fqcs.default <- function(x, argvals = NULL, rangeval  = NULL,...)
{

  m <- dim(x)[2] # quality characteristics

  if (is.null(argvals)) argvals <- seq(from = 1 ,to = m ,length.out = m )
  if (is.null(rangeval)) rangeval <- c(min(argvals),max(argvals))

  obj<-lab.fqcd(x = x, argvals = argvals, rangeval  = rangeval)

  result<-lab.fqcs.lab.fqcd(x = obj)

  return(result)
} #lab.fqcs

##' @rdname  lab.fqcs
##' @method lab.fqcs lab.fqcd
##' @inheritParams lab.fqcs.default
##' @export
lab.fqcs.lab.fqcd <- function(x, ...)
  #.........................................................................
  {
  if(is.null(x) || !inherits(x, "lab.fqcd"))
    stop("x must be an objects of class (or extending) 'lab.fqcd'")

  p <- x$p
  n <- x$n
  m <- x$m

  if(p ==1)
    stop("The object 'lab.fqcd' must be of length more than one laboratory")

  argvals <- x$curves.fdata[[1]]$argvals
  mean.i <- lapply(x$curves.fdata,func.mean)
  mean.i.m <- matrix(,nrow =p, ncol = m)
  for (i in 1:p)    mean.i.m[i,] <- mean.i[[i]]$data
  mean.i.lab.fqcd <- lab.fqcd(mean.i.m, argvals = argvals)
  
  s2.i <- lapply(x$curves.fdata,func.var)
  s2.i.m <- matrix(,nrow =p, ncol = m)
  for (i in 1:p)    s2.i.m[i,] <- s2.i[[i]]$data
  s2.i.lab.fqcd <- lab.fqcd(s2.i.m,argvals = argvals)
  
  curves <- matrix(,nrow = n*p, ncol = m)

  k <- 0
  for (i in 1:p){
    for (j in 1:n){
    k <- k+1
      curves[k,] <- x$curves[j,,i]
    }
  }
  fdata(curves, argvals = argvals)
  mean <- func.mean(fdata(curves, argvals = argvals))
  attr(mean, "object.name") <-"fdata" 
  S2 <- func.var(fdata(curves, argvals = argvals))
  attr(S2, "object.name") <-"fdata" 

result <- list (lab.fqcd = x, mean.i = mean.i.lab.fqcd, s2.i = s2.i.lab.fqcd,
                mean = mean, S2 = S2, p = p, n =n, m =m )

  oldClass(result) <- c("lab.fqcs")
  attr(result, "object.name") <- "lab.fqcs"
  attr(result, "type.data") <- "Functional Data List"
  
  return(result)
} # lab.fqcs
#.........................................................................
##' @rdname lab.fqcs
##' @method print lab.fqcs
##' @param x A \code{lab.fqcs} object for which a print is desired.
##' @export
print.lab.fqcs <- function(x, ...) str(x,1)
#.........................................................................
##' @rdname lab.fqcs
##' @method summary lab.fqcs
##' @param object A \code{lab.fqcs} object for which a summary is desired.
##' @export
summary.lab.fqcs <- function(object, ...)
  #.........................................................................
{

  type.data <- attributes(object)$object.name
 
  p <- object$p
  cat("\nNumber of laboratories: ", p)
  n <- object$n
  cat("\nNumber of replicates: ", n)
    
  result <- switch(type.data, 
                   "lab.fqcs" =  {
                     cat("\nSummary of fdata list:\n")
                     print(summary(object$lab.fqcd))
                     cat("\nSummary of mean fdata list:\n")
                     print(summary(object$mean.i))
                     cat("\nSummary of var fdata list:\n")
                     print(summary(object$s2.i))
                   } ,
                   "h.fqcs" =   {
                     cat("\nSummary of h statistics:\n")
                     print(summary(object$hi))
                   },
                   "k.fqcs" ={
                     cat("\nSummary of k statistics:\n")
                     print(summary(object$ki))
                   },
                   "QuantileWalter" = {
                     cat("\nSummary of Quantile by method's Walter:\n")
                     print(str(object))
                   },
                   "QuantileDepth" = {
                     cat("\nSummary of Quantile by method's Depth:\n")
                     print(str(object))
                   })
    #.........................................................................
} # summary.lab.fqcs


#-----------------------------------------------------------------------------#
# Main function to create a 'lab.fqcs' object
#-----------------------------------------------------------------------------#
##' This function is used to compute the FDA Mandel's h statistic.
##'
##' It develops an object of 'h.fqcs' class to perform statistical quality control analysis. 
##' This function is used to compute the functional approach of  Mandel's h statistic. 
##' It is specifically designed to deal with experimental data results defined by curves such as thermograms and spectra.

## @rdname  lab.fqcs
## @param x  Object lab.fqcd (Functional Quality Control Data)
##' @export
##' @references
##' \describe{
##'   \item{}{Febrero-Bande, M. and Oviedo, M. (2012),
##'    "Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28.}
##'   \item{}{Cuevas A., Febrero-Bande, M. and Fraiman, R. (2006), "On the use of the bootstrap for estimating functions with functional data". 
##'   Computational Statistics & Data Analysis 51, 2, 1063-1074. }
##'   \item{}{Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014), 
##'   "Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.}
##' }
##' @examples
##' library(ILS)
##' data(TG)
##' delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
##' curves.fqcd <- lab.fqcd(TG, argvals = delta)
##' fh <- h.fqcs(curves.fqcd)
##' xlab <- "Temperature(Grade Celsius)"
##' ylab <- "Mass (%)"
##' main <- "Functional hi  Estimation"
##' p <- fh$p
##' plot(fh,main = main, xlab = xlab, ylab = ylab,col=1:p,ylim=c(-3,3))
##' legend(10,3,paste("Lab",1:p),
##'       col=1:p,lty=1,lwd=2,cex=0.7)
h.fqcs <- function(x, ...) {
  UseMethod("h.fqcs")
}
##' @rdname h.fqcs
##' @method h.fqcs default
##' @inheritParams lab.fqcd
##' @param ... Arguments passed to or from methods.
h.fqcs.default <- function(x, argvals = NULL, rangeval  = NULL, ...)
{
  if (!is.array(x) & !is.data.frame(x) & !is.matrix(x) & !class(x)=="lab.fqcd")
    stop("object must be a array or a matrix or data.frame or lab.fqcd")
  
  m <- dim(x)[2] # quality characteristics

  if (is.null(argvals)) argvals <- seq(from = 1 ,to = m ,length.out = m )
  if (is.null(rangeval)) rangeval <- c(min(argvals),max(argvals))

  obj<-lab.fqcd(x = x, argvals = argvals, rangeval  = rangeval)

  result<-h.fqcs.lab.fqcd(x = obj)

  return(result)
} #h.fqcs

##' @rdname  h.fqcs
##' @method h.fqcs lab.fqcd
##' @inheritParams h.fqcs.default
##' @export


h.fqcs.lab.fqcd <- function(x, ...){


  if (!is.array(x) & !is.data.frame(x) & !is.matrix(x) & !class(x)=="lab.fqcd")
    stop("object must be a array or a matrix or data.frame or lab.fqcd")
  
  
  x.lab.fqcs <- lab.fqcs(x)
  p <- x$p
  n <- x$n
  m <- x$m
  argvals <- x$curves.fdata[[1]]$argvals
  mean.i <- x.lab.fqcs$mean.i[[1]]

  mean <- apply(mean.i, 2,mean)
  Sbar <- apply(mean.i, 2,sd)

  hs <- t(apply(mean.i, 1,function(x) {(x-mean)/Sbar}))

  hi.lab.fqcd <- lab.fqcd(hs, argvals = argvals)

  h.mean <- func.mean( hs)
  attr(h.mean, "object.name") <-"fdata"
  Sbar <-fdata(Sbar, argvals = argvals)
  attr(Sbar, "object.name") <-"fdata"
  result <- list(hi = hi.lab.fqcd, h.mean = h.mean, Sbar = Sbar, p = p, n =n, m =m)

  oldClass(result) <- c("lab.fqcs")
  attr(result, "object.name") <- "h.fqcs"
  attr(result, "type.data") <- "Functional Data List"
  return(result )
}

#-----------------------------------------------------------------------------#
# Main function to create a 'lab.fqcs' object
#-----------------------------------------------------------------------------#
##' This function is used to compute the FDA Mandel's k statistic 
##'
##' It develops an object of 'k.fqcs' class to perform statistical quality control analysis. 
##' This function is used to compute the functional approach of  Mandel's k statistic. 
##' It is specifically designed to deal with experimental data results defined by curves such as thermograms and spectra.

## @aliases k.fqcs summary.k.fqcs print.k.fqcs
## @rdname  lab.fqcs
## @param x  Object lab.fqcd (Functional Quality Control Data)
##' @export
##' @references
##' \describe{
##'   \item{}{Febrero-Bande, M. and Oviedo, M. (2012),
##'    "Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28.}
##'   \item{}{Cuevas A., Febrero-Bande, M. and Fraiman, R. (2006), "On the use of the bootstrap for estimating functions with functional data". 
##'   Computational Statistics & Data Analysis 51, 2, 1063-1074. }
##'   \item{}{Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014), 
##'   "Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.}
##' }
##' @examples
##' library(ILS)
##' data(TG)
##' delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
##' curves.fqcd <- lab.fqcd(TG, argvals = delta)
##' fk <- k.fqcs(curves.fqcd)
##' p <- fk$p
##' xlab <- "Temperature (Grade Celsius)"
##' ylab <- "Mass (%)"
##' main <- "Functional ki  Estimation"
##' plot(fk,main = main, xlab = xlab, ylab = ylab, ylim=c(0,3),col=1:p)
##' legend(10,3,paste("Lab",1:p),
##'       col=1:p,lty=1,lwd=2,cex=0.7)
k.fqcs <- function(x, ...) {
  UseMethod("k.fqcs")
}

##' @rdname k.fqcs
##' @method k.fqcs default
##' @inheritParams lab.fqcd
##' @param ... Arguments passed to or from methods.
k.fqcs.default <- function(x, argvals = NULL, rangeval  = NULL,...)
{
  if (!is.array(x) & !is.data.frame(x) & !is.matrix(x) & !class(x)=="lab.fqcd")
    stop("object must be a array or a matrix or data.frame or lab.fqcd")
  
  m <- dim(x)[2] # quality characteristics

  if (is.null(argvals)) argvals <- seq(from = 1 ,to = m ,length.out = m )
  if (is.null(rangeval)) rangeval <- c(min(argvals),max(argvals))

  obj<-lab.fqcd(x = x, argvals = argvals, rangeval  = rangeval)

  result<-k.fqcs.lab.fqcd(x = obj)

  return(result)
} #k.fqcs

##' @rdname  k.fqcs
##' @method k.fqcs lab.fqcd
##' @inheritParams k.fqcs.default
##' @export
k.fqcs.lab.fqcd <- function(x, ...){

  
  if (!is.array(x) & !is.data.frame(x) & !is.matrix(x) & !class(x)=="lab.fqcd")
    stop("object must be a array or a matrix or data.frame or lab.fqcd")
  
  x.lab.fqcs <- lab.fqcs(x)
  p <- x$p
  n <- x$n
  m <- x$m
  argvals <- x$curves.fdata[[1]]$argvals
  s2.i <- x.lab.fqcs$s2.i[[1]]
  s2r <- apply(s2.i, 2,mean)

  ks<- t(apply(s2.i, 1,function(x) {sqrt(x)/sqrt(s2r)}))

  ki.lab.fqcd <- lab.fqcd(ks, argvals = argvals)

  ks <-fdata(ks, argvals = argvals)
  k.mean <- func.mean( ks)
  attr(k.mean, "object.name") <- "fdata"
  s2r <-fdata(s2r, argvals = argvals)
  attr(s2r, "object.name") <- "fdata"
  result <- list(ki = ki.lab.fqcd, k.mean = k.mean, s2r = s2r, p = p, n =n, m =m)
  
  oldClass(result) <- c("lab.fqcs")
  attr(result, "object.name") <- "k.fqcs"
  attr(result, "type.data") <- "Functional Data List"
  
  return(result )
}
#-----------------------------------------------------------------------------#
# Main function to create a 'lab.fqcs' object
#-----------------------------------------------------------------------------#
##' It provides quantile estimates from FDA point of view.
##' 
##' It  develops a 'lab.fqcs' object to estimate functional quatiles by Walter's method (2011) from a pointwise point of view.
##' @rdname QuantileWalter 
##' @param x  Object of type fdata
##' @param quantile Probability defined in the interval [0,1]
##' @param central Logical argument. If  FALSE, functional quantile q is computed. 
##' If TRUE, two functional quantiles are obtained, those corresponding to curves (1-q / 2) and q / 2 .
##' @export
##' @references
##' \describe{
##'   \item{}{Febrero-Bande, M. and Oviedo, M. (2012),
##'    "Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28.}
##'   \item{}{Walter, S. (2011), Defining Quantiles for Functional Data: with an Application to the Reversal of Stock Price Decreases, 
##'   Department of Math. and Stat. The Uni. of Melbourne.}
##' }
##' @examples
##' \dontrun{
##' library(ILS)
##' data(TG)
##' delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
##' curves.fqcd <- lab.fqcd(TG, argvals = delta)
##' n <- curves.fqcd$n
##' m <- curves.fqcd$m
##' p <- curves.fqcd$p
##' curves.all <- TG[,,1]
##' for(i in 2:p) curves.all <- rbind(curves.all,TG[,,i])
##' curves.fdata <- fdata(mdata = curves.all,delta)
##' qw <- QuantileWalter(curves.fdata)
##' windows(20,10)
##' par(mfrow=c(1,2))
##' plot(qw, main="Quantiles of TG curves (95%)",col=c("red","blue"),lwd=2,legend = FALSE)
##' legend(50,80,c("Quantile 2.5%","Quantile 97.5%"),
##'       col=c("red","blue"),lty=c(1,1),lwd=1,cex=0.7)
##' plot(curves.fdata,main="Quantiles of TG curves (95%)",col="gray")
##' for(i in 1:2)
##' lines(qw[[i]],col="red",lty = 2,lwd = 2)
##' legend(50,80,c("Quantiles","TG Curves (105)"),
##'       col=c("red","gray"),lty=c(1,2),lwd=2,cex=0.7)
##' par(mfrow=c(1,1))
##' }
QuantileWalter <- function (x, quantile = 0.95, central = TRUE)
{
  if (!is.fdata(x))
    stop("object must be a fdata")

  m <- ncol(x[["data"]])
  x.low <- matrix(,ncol = m, nrow = 1)
  y.up <- matrix(,ncol = m, nrow = 1)
  argvals <- x$argvals

  if (central){

    x.low <- matrix(apply(x[["data"]], 2,
                                       quantile, probs = (1-quantile)/2, na.rm = TRUE), nrow = 1)

    x.low <-fdata(x.low, argvals = argvals)
    x.low$names$main <- "quantile.low "

    y.up <- matrix(apply(x[["data"]], 2,
                          quantile, probs = quantile+(1-quantile)/2, na.rm = TRUE), nrow = 1)

    y.up <-fdata(y.up, argvals = argvals)
    y.up$names$main <- "quantile.up "

    result <- list(low = x.low, up = y.up)

  }else{

    y.up <- matrix(apply(x[["data"]], 2,
                         quantile, probs = quantile, na.rm = TRUE), nrow = 1)

    y.up <-fdata(y.up, argvals = argvals)
    y.up$names$main <- "quantile"
    result <- list(up = y.up)
  }
  oldClass(result) <- c("lab.fqcs")
  attr(result, "object.name") <- "QuantileWalter"
  attr(result, "type.data") <- "Functional Data List"
  
  return (result)
}

#-----------------------------------------------------------------------------#
# Main function to create a 'lab.fqcs' object
#-----------------------------------------------------------------------------#
##' Creates a 'lab.fqcs' object to estimate functional quantiles using data depth procedures.
##' 
##' It defines a \code{lab.fqcs} object to estimate functional quantiles using data depth procedures (Lopez-Pintado and Romo, 2009).
##' The required functional quantiles are obtained from data depths values of each curve. 
##' If quantile argument is 0.9, 0.95 and 0.05 functional quantiles are obtained.
##' @rdname QuantileDepth
##' @param x  Object of type fdata
##' @param quantile Probability defined in the interval [0,1]
##' @export
##' @references
##' \describe{
##'   \item{}{Febrero-Bande, M. and Oviedo, M. (2012),
##'    "Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28.}
##'   \item{}{Lopez-Pintado, S. and Romo, J. (2009), "On the concept of depth for functional data", 
##'   Journal of the American Statistical Association, 104, 486-503. }
##' }
##' @examples
##' \dontrun{
##' library(ILS)
##' data(TG)
##' delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
##' curves.fqcd <- lab.fqcd(TG, argvals = delta)
##' n <- curves.fqcd$n
##' m <- curves.fqcd$m
##' p <- curves.fqcd$p
##' curves.all <- TG[,,1]
##' for(i in 2:p) curves.all <- rbind(curves.all,TG[,,i])
##' curves.fdata <- fdata(mdata = curves.all,delta)
##' qd <- QuantileDepth(curves.fdata)
##' windows(20,10)
##' par(mfrow=c(1,2))
##' plot(qd, main="Quantiles of TG curves (95%)",col=c("red","blue"),lwd=2,legend = FALSE)
##' legend(50,80,c("Quantile 2.5%","Quantile 97.5%"),
##'       col=c("red","blue"),lty=c(1,1),lwd=1,cex=0.7)
##' plot(curves.fdata,main="Quantiles of TG curves (95%)",col="gray")
##' for(i in 1:2)
##' lines(qd[[i]],col="red",lty = 2,lwd = 2)
##' legend(50,80,c("Quantiles","TG Curves (105)"),
##'       col=c("red","gray"),lty=c(1,2),lwd=2,cex=0.7)
##' par(mfrow=c(1,1))
##' }
QuantileDepth <- function (x, quantile = 0.95 )
{
  if (!is.fdata(x))
    stop("object must be a fdata")

  n <- nrow(x[["data"]])
  m <- ncol(x[["data"]])

  x.low <- matrix(,ncol = m, nrow = 1)
  y.up <- matrix(,ncol = m, nrow = 1)
  argvals <- x$argvals


  depth <- unlist(MBD(x = x[["data"]], plotting = FALSE)$MBD)
  index <- order(depth, decreasing=T )
  r <- ceiling( n*quantile)
  center <- x[["data"]][index[1:r], ]

  x.low <- matrix(apply(center, 2, min), nrow = 1)

  x.low <-fdata(x.low, argvals = argvals)
  x.low$names$main <- "quantile.low "

  y.up <- matrix(apply(center, 2, max), nrow = 1)

  y.up <-fdata(y.up, argvals = argvals)
  y.up$names$main <- "quantile.up "

  result <- list(low = x.low, up = y.up)
  oldClass(result) <- c("lab.fqcs")
  attr(result, "object.name") <- "QuantileDepth"
  attr(result, "type.data") <- "Functional Data List"
  
  return (result)
}

#-----------------------------------------------------------------------------#
# Main function to create a 'bootstrap.quantile' object
#-----------------------------------------------------------------------------#
##' Compute functional (FDA) Mandel's h and k statistics
##'
##' This function is used to compute functional (FDA)Mandel's h and k, statistics, 
##' required to perform Interlaboratory studies, and to detect non-consistent laboratories where data show a functional form (curve). In addition, 
##' bootstrap resampling methodology is used to estimate functional distributions. 
##' This allow to perform bootstrap confidence bands for FDA h and k statistics.
##'  
##' @param statistic Sample statistic used for the interlaboratory analysis.
##' By default, it uses sample h.
##' @param method Quantile method used to estimate the critical quantile of the h and k statistics.
##' @param quantile Probability with value in [0,1]
##' @param ball Logical argument. If draw = TRUE and ball = TRUE, i bootstrap curves and quantiles functions are plotted. 
##' They correspond to (1-alpha/2)*100 [\%]  most central  bootstrap resampling curves of q quantile. If draw = TRUE and ball = FALSE, 
##' the functional quantile q [\%] is determined.
##' @param alpha Significance level.
##' @param nb Number of bootstrap resamples.
##' @param smo Smoothing parameter for the bootstrap resamples, defined as a proportion of the sample variance matrix.
##' @param draw Default TRUE, it plots the bootstrap samples and the h or k statistic.  It depends on the ball parameter.
##' @param draw.control List that specifies the col, lty and lwd plot arguments for the objects lab.fqcs, statistic, IN and OUT.
##' @param x.co It speficies the x co-ordinates to be used to place a legend.
##' @param y.co It specifies the y co-ordinates to be used to place a legend.
##' @param legend Logical argument. Default is TRUE then The legend default is used.
##' @param col Color specifications
##' @param ... Arguments passed to or from methods.
##' @export
##' @references
##' \describe{
##'   \item{}{Febrero-Bande, M. and Oviedo, M. (2012),
##'    "Statistical computing in functional data analysis: the R package fda.usc". Journal of Statistical Software 51 (4), 1-28.}
##'   \item{}{Cuevas A., Febrero-Bande, M. and Fraiman, R. (2006), "On the use of the bootstrap for estimating functions with functional data". 
##'   Computational Statistics & Data Analysis 51, 2, 1063-1074. }
##'   \item{}{Naya, S., Tarrio-Saavedra. J., Lopez- Beceiro, J., Francisco Fernandez, M., Flores, M. and  Artiaga, R. (2014), 
##'   "Statistical functional approach for interlaboratory studies with thermal data". Journal of Thermal Analysis and Calorimetry, 118,1229-1243.}
##'   \item{}{Lopez-Pintado, S. and Romo, J. (2009), "On the concept of depth for functional data", 
##'   Journal of the American Statistical Association, 104, 486-503. }
##'   \item{}{Walter, S. (2011), Defining Quantiles for Functional Data: with an Application to the Reversal of Stock Price Decreases, 
##'   Department of Math. and Stat. The Uni. of Melbourne.}
##' }
##' @examples
##' \dontrun{
##' library(ILS)
##' data(TG)
##' delta <- seq(from = 40 ,to = 850 ,length.out = 1000 )
##' curves.fqcd <- lab.fqcd(TG, argvals = delta,rangeval = c(40,80))
##' draw.control = list(col = c("blue","grey"), 
##'                    lty = c(1, 1), lwd = c(2, 1))
##' #by Walter method
##' windows(20,10)
##' par(mfrow=c(1,2))
##' quantile95.w <- bootstrap.quantile(curves.fqcd, statistic = "h", 
##'                                   method = "Walter", smo = 0, 
##'                                   nb= 500, alpha = 0.05, quantile = 0.95,draw = TRUE,
##'                                   draw.control = draw.control,ylim=c(-3,3),x.co=50,y.co=3,
##'                                   main="Statistical h by the method's Walter")
##' }
bootstrap.quantile <- function(x, ...) {
  UseMethod("bootstrap.quantile")
}

##' @rdname bootstrap.quantile
##' @method bootstrap.quantile default
##' @inheritParams lab.fqcd
bootstrap.quantile.default <- function(x, argvals = NULL, rangeval  = NULL,
                           statistic = c("h","k"), method = c("Walter","Depth"),
                           alpha = 0.05, quantile = 0.9, ball = FALSE, nb = 200,
                           smo = 0, draw = TRUE, draw.control = NULL, x.co = NULL, y.co = NULL, 
                           legend = TRUE, col = NULL, ...)
{
  if (!is.array(x) & !class(x)=="lab.fqcd")
    stop("object must be a array or a matrix or data.frame or lab.fqcd")
  
  m <- dim(x)[2] # quality characteristics
  statistic <- match.arg(statistic)
  method <- match.arg(method)

  if (is.null(argvals)) argvals <- seq(from = 1 ,to = m ,length.out = m )
  if (is.null(rangeval)) rangeval <- c(min(argvals),max(argvals))

  obj<-lab.fqcd(x = x, argvals = argvals, rangeval  = rangeval)

  result<-bootstrap.quantile.lab.fqcd(x = obj,statistic = statistic, method = method,
                                           alpha = alpha, quantile = quantile, ball = ball, nb = nb,
                                           smo = smo, draw = draw, draw.control = draw.control, x.co = x.co, y.co = y.co, 
                                           legend = legend, col = col)

  return(result)
} #bootstrap.quantile

##' @rdname bootstrap.quantile
##' @method bootstrap.quantile lab.fqcd
##' @inheritParams bootstrap.quantile.default
##' @export
bootstrap.quantile.lab.fqcd <- function (x, statistic = c("h","k"), method = c("Walter","Depth"),
                                     alpha = 0.05, quantile = 0.9, ball = FALSE, nb = 200,
                            smo = 0, draw = TRUE, draw.control = NULL, x.co = NULL, y.co = NULL, 
                            legend = TRUE, col = NULL, ...)
{

  if (!is.array(x) & !class(x)=="lab.fqcd")
    stop("object must be a array or a matrix or data.frame or lab.fqcd")
  
  statistic <- match.arg(statistic)
  method <- match.arg(method)
  x<- lab.fqcs(x)
  p <- length(x$lab.fqcd$curves.fdata)
  data <- x$lab.fqcd$curves.fdata[[1]]$data
  nr <- nrow(data)
  nc <- ncol(data)
  tt <- x$lab.fqcd$curves.fdata[[1]]$argvals
  rtt <- x$lab.fqcd$curves.fdata[[1]]$rangeval
  names <- x$lab.fqcd$curves.fdata[[1]]$names

  if (is.null(x.co)) x.co <- min(tt)
  if (is.null(y.co)) y.co <- 0.99 * max(data)
  if (is.null(col))  col <- terrain.colors(p)

  m.data <- array(,dim = c(nr,nc,p))
  curves <- matrix(,nrow = nr*p, ncol = nc)

  k <- 0
  for (i in 1:p){
    for (j in 1:nr){
      k <- k+1
      curves[k,] <- x$lab.fqcd$curves[j,,i]
    }
  }
# statistic sample
  if (statistic =="h"){
    m.h <- h.fqcs(x$lab.fqcd)
    m.his <- matrix(,nrow = p, ncol = nc)
    for (l in 1:p){
      m.his[l,] <- as.vector(m.h$hi$curves.fdata[[l]]$data)
    }
    m.his <-fdata(m.his, argvals = tt)
    data <- m.his$data
    if (method == "Walter"){
      m.quantile <- QuantileWalter(m.his, quantile = quantile, central = TRUE)
    }else{
      m.quantile <- QuantileDepth(m.his, quantile = quantile)
    }
    m.statistic <- m.h$hi

  }else{
    m.k <- k.fqcs(x$lab.fqcd)
    m.kis <- matrix(,nrow = p, ncol = nc)
    for (l in 1:p){
      m.kis[l,] <- as.vector(m.k$ki$curves.fdata[[l]]$data)
    }
    m.kis <-fdata(m.kis, argvals = tt)
    data <- m.kis$data
      m.quantile <- QuantileWalter(m.kis, quantile = quantile, central = FALSE)
    m.statistic <- m.k$ki
  }

  pb = txtProgressBar(min = 0, max = nb, style = 3)
  distboot <- matrix(NA, nrow = nb)
  estboot.low <- matrix(NA, nrow = nb, ncol = nc)
  estboot.up <- matrix(NA, nrow = nb, ncol = nc)

for (i in 1:nb) {
  setTxtProgressBar(pb, i - 0.5)
    for (j in 1:p ){
      m.data[,,j] <- curves[sample(1:nr, size = nr, replace = TRUE),]
      if (smo > 0)
        m.data[,,j] <- m.data[,,j] + mvrnorm(n = nr, rep(0, nc), var(x[[1]][,,j]) * smo)
     }

    b.fqcd <- lab.fqcd(m.data, argvals = tt)

    if (statistic =="h"){
      boot.h <- h.fqcs(b.fqcd)
      b.his <- matrix(,nrow = p, ncol = nc)
      for (l in 1:p){
        b.his[l,] <- as.vector(boot.h$hi$curves.fdata[[l]]$data)
      }
      b.his <-fdata(b.his, argvals = tt)
      if (method == "Walter"){
        boot.quantile <- QuantileWalter(b.his, quantile = quantile, central = TRUE)
      }else{
        boot.quantile <- QuantileDepth(b.his, quantile = quantile)
      }
      estboot.low[i, ] <- boot.quantile$low$data
      estboot.up[i, ] <- boot.quantile$up$data
    }
    else{
      boot.k <- k.fqcs(b.fqcd)
      b.kis <- matrix(,nrow = p, ncol = nc)

      for (l in 1:p){
        b.kis[l,] <- as.vector(boot.k$ki$curves.fdata[[l]]$data)
      }

      b.kis <-fdata(b.kis, argvals = tt)
        boot.quantile <- QuantileWalter(b.kis, quantile = quantile, central = FALSE)

      estboot.up[i, ] <- boot.quantile$up$data

    }

  setTxtProgressBar(pb, i)
}

  close(pb)


  if (statistic == "h"){
    center.low <- func.mean(fdata(estboot.low,tt))
    distboot.low <- metric.lp(fdata(estboot.low, tt, rtt), center.low)
    dist.low <- max(distboot.low[rank(distboot.low) <= floor((1 - alpha)*nb)])
    resample.low <- fdata(estboot.low, tt, rtt, names)

    center.up <- func.mean(fdata(estboot.up,tt))
    distboot.up <- metric.lp(fdata(estboot.up, tt, rtt), center.up)
    dist.up <- max(distboot.up[rank(distboot.up) <= floor((1 - alpha)*nb)])
    resample.up <- fdata(estboot.up, tt, rtt, names)
  }
  else{

    center.up <- func.mean(fdata(estboot.up,tt))
    distboot.up <- metric.lp(fdata(estboot.up, tt, rtt), center.up)
    dist.up <- max(distboot.up[rank(distboot.up) <= floor((1 - alpha)*nb)])
    resample.up <- fdata(estboot.up, tt, rtt, names)

  }

if (draw) {
  if (is.null(draw.control))
    draw.control = list(col = c("blue","grey" ),
                        lty = c(1, 1), lwd = c(1, 2))
  if (is.null(draw.control$lwd))
    draw.control$lwd = c(1, 2)
  if (is.null(draw.control$lty))
    draw.control$lty = c(1, 1)
  if (is.null(draw.control$col))
    draw.control$col = c("blue","grey" )

  plot.lab.fqcd(m.statistic, legend = FALSE, ...)


  if (statistic == "h"){

    if (ball){

      lines(center.low, lwd = draw.control$lwd[1], lty = draw.control$lty[1],
            col = draw.control$col[1])
      lines(center.up, lwd = draw.control$lwd[1], lty = draw.control$lty[1],
            col = draw.control$col[1])


      lines(resample.low[distboot.low <= dist.low], lwd = draw.control$lwd[2],
              lty = draw.control$lty[2], col = draw.control$col[2])

      lines(resample.up[distboot.up <= dist.up], lwd = draw.control$lwd[2],
              lty = draw.control$lty[2], col = draw.control$col[2])

      if (legend ==TRUE){
        legend(x = x.co, y = y.co, legend = c(paste("Lab",1:p),
                                              "Quantiles", "Bootstrap curves IN"),
               lty = draw.control$lty,
               lwd = draw.control$lwd, col = c(col,draw.control$col), cex = 0.9,
               box.col = 1)
      }
    }else{

      lines(center.low, lwd = draw.control$lwd[1], lty = draw.control$lty[1],
            col = draw.control$col[1])
      lines(center.up, lwd = draw.control$lwd[1], lty = draw.control$lty[1],
            col = draw.control$col[1])

      if (legend ==TRUE){
      legend(x = x.co, y = y.co, legend = c(paste("Lab",1:p),
                                            "Quantiles"),
             lty = draw.control$lty[1],
             lwd = draw.control$lwd[1], col = c(col,draw.control$col[1]), cex = 0.9,
             box.col = 1)
      }
    }

  }else{

      if (ball){

        lines(resample.up[distboot.up <= dist.up], lwd = draw.control$lwd[2],
              lty = draw.control$lty[2], col = draw.control$col[2])

        lines(center.up, lwd = draw.control$lwd[1], lty = draw.control$lty[1],
              col = draw.control$col[1])

        if (legend ==TRUE){
        legend(x = x.co, y = y.co, legend = c(paste("Lab",1:p),
                                              "Quantiles", "Bootstrap curves IN"),
               lty = draw.control$lty,
               lwd = draw.control$lwd, col = c(col,draw.control$col), cex = 0.9,
               box.col = 1)
        }

      }else{
        lines(center.up, lwd = draw.control$lwd[1], lty = draw.control$lty[1],
              col = draw.control$col[1])

        
        if (legend ==TRUE){
        legend(x = x.co, y = y.co, legend = c(paste("Lab",1:p),
                                              "Quantiles"),
               lty = draw.control$lty[1],
               lwd = draw.control$lwd[1], col = c(col,draw.control$col[1]), cex = 0.9,
               box.col = 1)
        }

      }


  }

}
  
if  (statistic == "h"){ 
result <- list(quantile.low = center.low, dband.low = dist.low, rep.dist.low = distboot.low,
            resample.low = resample.low,
                                    m.statistic = m.statistic,
            quantile.up = center.up, dband.up = dist.up, rep.dist.up = distboot.up,
            resample.up = resample.up, p = p, n = nr)
}else{
  result <- list(m.statistic = m.statistic,
                 quantile.up = center.up, dband.up = dist.up, rep.dist.up = distboot.up,
                 resample.up = resample.up, p = p, n = nr, statistic = statistic)
}
 oldClass(result) <- c("bootstrap.quantile")
 attr(result, "object.name") <- "QuantileBootstrap"
 attr(result, "type.data") <- "Functional Data List"
 
 return(result)
}

##' @rdname  bootstrap.quantile
##' @method print bootstrap.quantile
##' @param x A \code{bootstrap.quantile} object for which a print is desired.
##' @export
print.bootstrap.quantile <- function(x, ...) str(x,1)
#.........................................................................

##' @rdname  bootstrap.quantile
##' @method summary bootstrap.quantile
##' @param object A \code{bootstrap.quantile} object for which a summary is desired.
##' @export

summary.bootstrap.quantile <- function(object, ...)
  #.........................................................................
{
  p <- object$p
  cat("\nNumber of laboratories: ", p)
  
  n <- object$n
  cat("\nNumber of replicates: ", n)
  
  cat("\nSummary of statistical sampling:\n")
  print(summary(object["m.statistic"]))
  
  if (object["statistic"]=="h"){
  cat("\nSummary of quantile low:\n")
  print(summary(object["quantile.low"]))
  
  cat("\nSummary of quantile up:\n")
  print(summary(object["quantile.up"]))
  }else{
    cat("\nSummary of quantile:\n")
    print(summary(object["quantile.up"]))
  }
  
  #.........................................................................
} # summary.bootstrap.quantile
