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
# Main function to create a 'ils.fqcs' object
#-----------------------------------------------------------------------------#
##' It developes an object of class 'ils.fqcs'
##' 
##' Create an object of class 'ils.fqcs' to perform statistical quality control. 
##' This function is used to compute requested FDA.
##'
## @param x  Object ils.fqcdata (Functional Quality Control Data)
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
##' fqcdata <- ils.fqcdata(TG, p = 7, argvals = delta)
##' xlab <- "Temperature/ C"
##' ylab <- "Mass/ %"
##' fqcstat <- ils.fqcs(fqcdata)
##' plot(fqcstat, xlab = xlab, ylab = ylab,legend = TRUE) 

ils.fqcs <- function(x, ...) {
  UseMethod("ils.fqcs")
}

##' @rdname ils.fqcs
##' @method ils.fqcs default
##' @inheritParams ils.fqcdata
##' @param ... Arguments passed to or from methods.
ils.fqcs.default <- function(x, argvals = NULL, rangeval  = NULL,...)
{

  m <- dim(x)[2] # quality characteristics

  if (is.null(argvals)) argvals <- seq(from = 1 ,to = m ,length.out = m )
  if (is.null(rangeval)) rangeval <- c(min(argvals),max(argvals))

  obj<-ils.fqcdata(x = x, argvals = argvals, rangeval  = rangeval)

  result<-ils.fqcs.ils.fqcdata(x = obj)

  return(result)
} #ils.fqcs

##' @rdname  ils.fqcs
##' @method ils.fqcs ils.fqcdata
##' @inheritParams ils.fqcs.default
##' @export
ils.fqcs.ils.fqcdata <- function(x, ...)
  #.........................................................................
  {
 
  if(is.null(x) || !inherits(x, "ils.fqcdata"))
    stop("x must be an objects of class (or extending) 'ils.fqcdata'")
  
  
  p <- x$p
  n <- x$n
  m <- x$m
  index <- rep(x$index.laboratory,each = n)
  
  fdata.m <- x$ils.fdata
  S2 <- func.var(fdata.m)
  mean <- func.mean(fdata.m)
  fac <- factor(index)
  
  ldata <- list("df"=data.frame(fac),"fdataobj"=fdata.m)
  
  mean.i <- func.ils.formula(fdataobj~fac, data = ldata, func = func.mean)
  S2.i <- func.ils.formula(fdataobj~fac, data = ldata, func = func.var)
  Sbar <- fdata(apply(mean.i[["data"]], 2,sd))
  S2r <- fdata(apply(S2.i[["data"]], 2,mean))
  
  #Statistics
  ils.h <- fdata(t(apply(mean.i[["data"]], 1,function(x) 
  {(x-mean[["data"]])/Sbar[["data"]]})))
  ils.k<- fdata(t(apply(S2.i[["data"]], 1,function(x)
  {sqrt(x)/sqrt(S2r[["data"]])})))
  
  d_h.m <- c(norm.fdata(ils.h))
  names(d_h.m) <- paste("lab",1:p)
  
  d_k.m <- c(norm.fdata(ils.k))
  names(d_k.m) <- paste("lab",1:p)   
  
  result <- list(fqcdata = x, ils.h = ils.h, ils.k = ils.k, d_h.m = d_h.m, d_k.m = d_k.m, mean = mean, S2 = S2, 
                 mean.i = mean.i, S2.i = S2.i, Sbar = Sbar, S2r = S2r, 
                 p = p, n = n, m = m)

  oldClass(result) <- c("ils.fqcs")

  return(result)
} # ils.fqcs

#-----------------------------------------------------------------------------#
##' Descriptive measures for functional data.
##' 
##' Central and dispersion measures for functional data.
#.........................................................................
##' @rdname func.ils.formula
##' @param formula a formula, such as y ~ group, where y is a fdata object 
##' to be split into groups according to the grouping variable group (usually a factor).
##' @param data List that containing the variables in the formula. The item called "df" 
##' is a data frame with the grouping variable. The item called "y" is a fdata object.
##' @param drop logical indicating if levels that do not occur should be dropped (if f is a factor or a list).
##' @param func  Measures for functional data.
##' @export
#.........................................................................
func.ils.formula <- function (formula, data = NULL, drop = FALSE, func = func.mean) 
{
  tf <- terms.formula(formula)
  fac <- attr(tf, "term.labels")
  if (attr(tf, "response") > 0) 
    response <- as.character(attr(tf, "variables")[2])
  if (missing(formula) || (length(formula) != 3L)) 
    stop("'formula' missing or incorrect")
  ldata <- data
  data <- ldata$df
  if (is.null(ldata$df)) 
    stop("'df' element is missing or incorrect")
  if (is.vector(data)) 
    data <- as.data.frame(data)
  if (is.matrix(data)) 
    data <- as.data.frame(data)
  f <- ldata$df[[fac]]
  dat <- ldata[[response]]
  if (!is.factor(f)) 
    f <- factor(f)
  nlev <- nlevels(f)
  lev <- levels(f)
  if (is.matrix(dat$data)) 
    dat$data <- data.frame(dat$data)
  out <- split(dat$data, f, drop = drop)
  out2 <- func(fdata(out[[lev[1]]], dat$argvals, dat$rangeval, 
                     dat$names))
  for (i in 2:nlev) out2 <- c(out2, func(fdata(out[[lev[i]]], 
                                               dat$argvals, dat$rangeval, dat$names)))
  rownames(out2$data) <- lev
  return(out2)
}
#.........................................................................
#.........................................................................
##' @rdname ils.fqcs
##' @method print ils.fqcs
##' @param x A \code{ils.fqcs} object for which a print is desired.
##' @export
print.ils.fqcs <- function(x, ...) str(x,1)
#.........................................................................
##' @rdname ils.fqcs
##' @method summary ils.fqcs
##' @param object A \code{ils.fqcs} object for which a summary is desired.
##' @export
summary.ils.fqcs <- function(object, ...)
  #.........................................................................
{

  if(is.null(object) || !inherits(object, "ils.fqcs"))
    stop("x must be an objects of class (or extending) 'ils.fqcs'")
  
  type.data <- class(object)
 
  p <- object$p
  cat("\nNumber of laboratories: ", p)
  n <- object$n
  cat("\nNumber of replicates: ", n)
    
  result <- switch(type.data, 
                   "lab.fqcs" =  {
                     cat("\nSummary of lab.fqcd:\n")
                     print(summary(object$lab.fqcd))
                     cat("\nSummary of mean fdata:\n")
                     print(summary(list(mean.i = object$mean.i)))
                     cat("\nSummary of variance fdata:\n")
                     print(summary(list(s2.i = object$s2.i)))
                   } ,
                   "h.fqcs" =   {
                     cat("\nSummary of h statistics:\n")
                     print(summary(object$hi))
                   },
                   "k.fqcs" ={
                     cat("\nSummary of k statistics:\n")
                     print(summary(object$ki))
                   })
    #.........................................................................
} # summary.ils.fqcs


 #-----------------------------------------------------------------------------#
 # Main function to create a 'mandel.fqcs' object
 #-----------------------------------------------------------------------------#
 ##' This function is used to compute the FDA Mandel's h and k statistic.
 ##'
 ##' It develops an object of 'mandel.fqcs' class to perform statistical quality control analysis.
 ##' This function is used to compute the functional approach of  Mandel's h and k statistic.
 ##' It is specifically designed to deal with experimental data results defined by curves such as thermograms and spectra.

 ##' @rdname  mandel.fqcs
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
 ##' fqcdata <- ils.fqcdata(TG, p = 7, argvals = delta)
 ##' mandel.tg <- mandel.fqcs(fqcdata.tg,nb = 200)
 ##' plot(mandel.tg,legend = F,col=c(rep(3,5),1,1))
 ##'   }
 mandel.fqcs <- function(x, ...) {
   UseMethod("mandel.fqcs")
 }

##' @rdname mandel.fqcs
##' @method mandel.fqcs default
##' @inheritParams ils.fqcdata
mandel.fqcs.default <- function(x, p = NULL, index.laboratory = NULL, argvals = NULL, rangeval = NULL,
                               names = NULL, ...)
{
 if (!is.array(x) & !is.data.frame(x) & !is.matrix(x) & !class(x)=="ils.fqcdata")
   stop("object must be a array or a matrix or data.frame or ils.fqcdata")



 obj<-ils.fqcdata(x = x, argvals = argvals, rangeval  = rangeval, names = names)

 result<-mandel.fqcs.ils.fqcdata(x = obj)

 return(result)
} #mandel.fqcs

##' @rdname  mandel.fqcs
##' @method mandel.fqcs ils.fqcdata
##' @param fdep Type of depth measure, by default depth.mode.
##' @param outlier = TRUE
##' @param trim The alpha of the trimming.
##' @param alpha Significance level, by defaul 1\%.
##' @param nb The number of bootstrap samples.
##' @param smo The smoothing parameter for the bootstrap samples.
##' @param ... Further arguments passed to or from other methods.
##' @export
mandel.fqcs.ils.fqcdata<- function(x, fdep = depth.mode, outlier = TRUE, trim = 0.01,
                                  alpha = 0.01, nb = 200, smo = 0.05, ...){


 if (!is.array(x) & !is.data.frame(x) & !is.matrix(x) & !class(x)=="ils.fqcdata")
   stop("object must be a array or a matrix or data.frame or ils.fqcdata")

 fdata.m <- x$ils.fdata
 p <- x$p
 n <- x$n
 index <- rep(x$index.laboratory,each = n)
 mdata <- x$m
 alpha = alpha/p

 fqcs <- ils.fqcs(x)
 n.h <- fqcs$d_h.m
 n.k <- fqcs$d_k.m

 if(outlier == TRUE){

   out <- outliers.ils(x = fdata.m, fdep = fdep,trim = trim)

   if (length(out$outliers)!=0){
     fdata.in <- fdata.m[-as.numeric(out[[1]])]
     fdata.out <- fdata.m[as.numeric(out[[1]])]
   }else{
     fdata.in <- fdata.m
     fdata.out <- NULL
   }
 }else{
   fdata.in <- fdata.m
   fdata.out <- NULL
 }

 boot.sim <- boot.sim.set(fdata.in, smo = smo)
 ils.h.b <- NULL
 ils.k.b <- NULL
 d_h.b <- matrix(NA,nrow = nb, ncol = p)
 d_k.b <- matrix(NA,nrow = nb, ncol = p )

 for (b in seq_len(nb)) {

   boot.samples <- ils.fqcdata(t(boot.sim(rep = n * p)),p = p)
   #plot(boot.samples)

   boot.fqcs <- ils.fqcs(boot.samples)
   # h
   ils.h.b <- rbind(ils.h.b , boot.fqcs$ils.h[["data"]])
   d_h.b[b,] <- boot.fqcs$d_h.m
   # k
   ils.k.b <- rbind(ils.k.b , boot.fqcs$ils.k[["data"]])
   d_k.b[b,] <- boot.fqcs$d_k.m

 } # for (b in seq_len(nb))

 # h
 distboot.h <- c(d_h.b)
 dist.h <- max(distboot.h[rank(distboot.h) <= floor((1 - alpha)* (nb*p))])

 ind.h <- which(distboot.h <= dist.h)
 fmax.h <- fdata(apply(ils.h.b[ind.h,],2,max))
 fmin.h <- fdata(apply(ils.h.b[ind.h,],2,min))

 # k
 distboot.k <- c(d_k.b)
 dist.k <- max(distboot.k[rank(distboot.k) <= floor((1 - alpha)* (nb*p))])

 ind.k <- which(distboot.k <= dist.k)
 fmax.k <- fdata(apply(ils.k.b[ind.k,],2,max))
 fmin.k <- fdata(apply(ils.k.b[ind.k,],2,min))


 result<- list(m.ils.fdata = fdata.m,
               m.fest.h = fqcs$ils.h, m.fest.k = fqcs$ils.k,
               norm.h = n.h, norm.k = n.k,
               fdata.in = fdata.in,
               fdata.out = fdata.out, distboot.h = distboot.h,
               fmax.h = fmax.h, fmin.h = fmin.h, dist.h = dist.h,
               distboot.k = distboot.k, fmax.k = fmax.k, fmin.k = fmin.k, dist.k = dist.k,
               ind.h = which(fqcs$d_h.m >= dist.h) , ind.k = which(fqcs$d_k.m >= dist.k), fqcs = fqcs)

 oldClass(result) <- c("mandel.fqcs")

 return(result )
}
#-----------------------------------------------------------------------------
## @rdname boot.sim.set
##' Bootstrap samples of a functional statistic
##'
##' \code{data.bootstrap} provides bootstrap samples for functional data.
##' @param x A \code{fdata} object.
##' @param smo The smoothing parameter for the bootstrap samples.
##' @export
boot.sim.set <- function(x, smo = 0.05){
 data <- x[["data"]]
 mdata <- ncol(data)
 if (smo > 0) {
   # L.cov <- t(chol(var(data)* smo))
   # A Choleski decomposition might be faster, but the eigendecomposition is stabler...
   eS <- eigen(var(data)* smo, symmetric = TRUE)
   L.cov <- with(eS, vectors %*% diag(sqrt(pmax(values, 0)), mdata))
 }

 boot.sim <- function(rep){
   boot.data <- data[sample(seq_len(nrow(data)), size = rep, replace = TRUE),]
   if (smo > 0) {
     err.norm <- matrix(rnorm(mdata * rep), nrow = mdata)
     data.err <- L.cov %*% err.norm  #simulacion de errores Y = L %*% Z (usando la matriz de covarianzas de Y)
     boot.data <- boot.data + t(data.err)
   }
   boot.data <- if (rep == 1) drop(boot.data)  else t(boot.data)
   return(boot.data)
 }
 return(boot.sim)
}
#.........................................................................

##' Detecting outliers for functional dataset
##'
##' Procedure for detecting funcitonal outliers.
##' @param x A \code{fdata}
##' @param fdep Type of depth measure, by default depth.mode.
##' @param trim The percentaje of the trimming, by default is 1\%
##' @export
outliers.ils <- function(x, fdep = depth.FM, trim = 0.01){
if (!is.numeric(trim) || length(trim) != 1L || !(trim >= 0 & trim < 1))
 stop("'trim' must be a numeric scalar between 0 and 1")

#fdep <- match.arg(fdep)
if (!is.fdata(x)) x = fdata(x)

n <- nrow(x)
m <- ncol(x)

if (is.null(n) && is.null(m))
 stop("ERROR IN THE DATA DIMENSIONS")

d <-  fdep(x)$dep
index <- order(d, decreasing = FALSE)
num <- floor(trim * n)
outliers <- index[1:num]
good <- index[min(num + 1, n) : n]
muestra.trim <- x[good, ]

res <- list(outliers = outliers, depth = d, muestra.trim = muestra.trim, n.out = num)

return(res)

}

