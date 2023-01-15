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
# Main function to create a 'lab.qcs' object
#-----------------------------------------------------------------------------#
##' Create an object of class 'lab.qcs' to perform statistical quality control.
##' This function is used to compute statistics required for plotting Statitics
##'
##' It develops an object of \code{lab.qcs}-code{link{class}} to perform statistical quality control.
##' This function is used to compute the requested statistics to be summarized and ploted.
##'
## @aliases lab.qcs summary.lab.qcs print.lab.qcs
##' @param x  An object of class \code{lab.qcdata} (Univariate Quality Control Data).
##' @param ... Other arguments passed to or from methods.
##' @export
##' @examples
##'
##' library(ILS)
##' data(Glucose)
##' Glucose.qcdata <- lab.qcdata(Glucose)
##' str(Glucose.qcdata)
##' Glucose.qcs <- lab.qcs(Glucose.qcdata)
##' str(Glucose.qcs)
##' summary(Glucose.qcs)
lab.qcs <- function(x, ...)
  #.........................................................................
  {
  if(is.null(x) || !inherits(x, "lab.qcdata"))
    stop("x must be an objects of class (or extending) 'lab.qcdata'")

  p <- length(unique(x$laboratory))
  m <- length(unique(x$material))
  n <- length(unique(x$replicate))
  material<-unique(x$material)
  laboratory<-unique(x$laboratory)

  stat.material <- data.frame(mean = vector(,length = m),
                           S = vector(,length = m),
                           S_r = vector(,length = m),
                           S_B = vector(,length = m),
                           S_R = vector(,length = m))

  statistics <- data.frame(laboratory = vector(,length = p*m),
                           material = vector(,length = p*m),
                           mean.i = vector(,length = p*m),
                           s.i = vector(,length = p*m))

  data <- x$x

    statistics[,1] <- as.factor(rep(laboratory,each = m))
    statistics[,2] <- as.factor(rep(material,p))
    statistics[,3] <- c(tapply(data,list(x$material,x$laboratory),mean))
    statistics[,4] <- c(tapply(data,list(x$material,x$laboratory),sd))

    stat.material[,1] <- tapply(statistics$mean.i,statistics$material,mean)
    stat.material[,2] <- tapply(statistics$s.i,statistics$material,sd)
    f.S_r <- function(s.i) {sqrt(mean(s.i^2))}
    S_r <- stat.material[,3] <- tapply(statistics$s.i,statistics$material,f.S_r)

    S_B <- stat.material[,4] <- tapply(statistics$mean.i,statistics$material,sd)
    stat.material[,5] <- sqrt(S_B^2 + ((n-1)/n)*S_r^2)

    rownames(stat.material) <- material

 result <- list (lab.qcdata = x, statistics.Laboratory = statistics,
                 statistics.material =  stat.material, p = p, n = n, m = m )
 oldClass(result)<-c("lab.qcs")
 attr(result, "object.name") <- attributes(x)$data.name
 attr(result, "type.data") <- "lab.qcs"

return(result)
} # lab.qcs
#.........................................................................

##' @rdname lab.qcs
##' @method print lab.qcs
##' @param x An object of class \code{lab.qcs} for which a print is desired.
##' @export
print.lab.qcs <- function(x, ...) str(x,1)
#.........................................................................
##' @rdname lab.qcs
##' @method summary lab.qcs
##' @param object An object of class \code{lab.qcs} for which a summary is desired.
##' @export
summary.lab.qcs <- function(object, ...)
  #.........................................................................
{


  type.data <- attributes(object)$type.data

  cat("\nNumber of laboratories: ", object$p)
  cat("\nNumber of materials: ", object$m)
  cat("\nNumber of replicate: ", object$n)

  result <- switch(type.data,
                  "lab.qcs" =  {
                    cat("\nSummary for Laboratory (means):\n")
                    st <- with(object$lab.qcdata,
                               tapply(x,
                                      list(material,
                                           laboratory), mean))
                    print(st)

                    cat("\nSummary for Laboratory (Deviations):\n")
                    st <- with(object$lab.qcdata,
                               tapply(x,
                                      list(material,
                                           laboratory), sd))
                    print(st)

                    cat("\nSummary for Material:\n")
                    print(object$statistics.material)

                    },
                  "h.qcs" =   {

                    cat("\nCritical value: ", object[[7]])
                    cat("\nBeyond limits of control:", "\n")
                    print(object[[8]])
                  },
                  "k.qcs" ={

                    cat("\nCritical value: ", object[[7]])
                    cat("\nBeyond limits of control:", "\n")
                    print(object[[8]])
                    })

  invisible()
  #.........................................................................
} # summary.qcs

#-------------------------------------------------------------------------
# Statistic h
#-------------------------------------------------------------------------
##' Function to estimate the univariate Mandel's  h statistic
##'
##' This function is used to compute the Mandel's h statistic.
##' @param x   An object of class \code{lab.qcdata} (Univariate Quality Control Data).
##' @export
##' @references
##' Wilrich Peter-T. (2013),  Critical values of Mandel's h and k,
##' the Grubbs and the Cochran test statistic. Asta-Advances in Statistical Analysis, 97(1):1-10.
##' 
##' ASTM E 691 (1999), Standard practice for conducting an interlaboratory study
##'   to determine the precision of a test method. American Society for Testing and Materials. West Conshohocken, PA, USA.
##' @examples
##'
##' library(ILS)
##' data(Glucose)
##' Glucose.qcdata <- lab.qcdata(Glucose)
##' str(Glucose.qcdata)
##' h<- h.qcs(Glucose.qcdata, alpha = 0.005)
##' summary(h)
##' plot(h)

h.qcs <- function(x, ...) {
  UseMethod("h.qcs")
}

##' @rdname h.qcs
##' @method h.qcs default
##' @inheritParams lab.qcdata
##' @param alpha The significance level (0.05 by default)
##' @param ... Other arguments passed to or from methods.
h.qcs.default <- function(x, var.index=1,replicate.index  =  2, material.index  =  3,
                          laboratory.index=4,  data.name = NULL, alpha = 0.05, ...)
  {
  if (is.null(data.name)) data.name <- "Statistical Mandel h"

    obj<-lab.qcdata(data = x, var.index=var.index,replicate.index  =  replicate.index,
                 material.index  =  material.index,
                 laboratory.index=laboratory.index,  data.name = data.name)

    result<-h.qcs.lab.qcdata(x = obj,  alpha = alpha)

  return(result)
} #h.qcs

##' @rdname  h.qcs
##' @method h.qcs lab.qcdata
##' @inheritParams h.qcs.default
##' @export
h.qcs.lab.qcdata <- function(x, alpha = 0.05, ...)
{

  if(is.null(x) || !inherits(x, "lab.qcdata"))
    stop("x must be an objects of class (or extending) 'lab.qcdata'")

  data.name <- attributes(x)$data.name
  x.lab.qcs <- lab.qcs(x)

  statistics <- x.lab.qcs$statistics.material
  mean.i <- x.lab.qcs$statistics.Laboratory$mean.i
  p <- x.lab.qcs$p
  n <- x.lab.qcs$n
  m <- x.lab.qcs$m

  hcrit <- (p-1)*qt((1-alpha/2),(p-2))/sqrt(p*(p-2+qt((1-alpha/2),(p-2))^2))

  material <- row.names(x.lab.qcs$statistics.material)
  laboratory <- unique(x.lab.qcs[[1]]$laboratory)
  h.i <- matrix(,nrow = p,ncol = m)
  for(i in 1:m)
  {
    ind <- x.lab.qcs$statistics.Laboratory$material==material[i]
    h.i[,i] <- (mean.i[ind]-statistics$mean[i])/statistics$S[i]
  }



  colnames(h.i) <- material
  rownames(h.i) <- laboratory
  violations <- abs(h.i) <= hcrit
  result <- list (lab.qcdata = x, lab.qcs = x.lab.qcs, p = p, n = n, m = m,
                  h = h.i, h.critial = hcrit, violations = violations, data.name = data.name )

  oldClass(result) <- c("lab.qcs")
  attr(result, "object.name") <- data.name
  attr(result, "type.data") <- "h.qcs"


  return(result)

} #h.qcs

# Statistic k
#-------------------------------------------------------------------------
##' Function to calcute the Mandel's k statistic
##'
##' This function is used to compute the statistic k of Mandel.
##' @param x   An object of class \code{lab.qcdata} (Univariate Quality Control Data).
##' @export
##' @references
##' Wilrich Peter-T. (2013),  Critical values of Mandel's h and k,
##' the Grubbs and the Cochran test statistic. Asta-Advances in Statistical Analysis, 97(1):1-10.
##' 
##' ASTM E 691 (1999), Standard practice for conducting an interlaboratory study
##' to determine the precision of a test method. American Society for Testing and Materials. West Conshohocken, PA, USA.
##' @examples
##'
##' library(ILS)
##' data(Glucose)
##' Glucose.qcdata <- lab.qcdata(Glucose)
##' str(Glucose.qcdata)
##' k<- k.qcs(Glucose.qcdata, alpha = 0.005)
##' summary(k)
##' plot(k)
k.qcs <- function(x, ...) {
  UseMethod("k.qcs")
}

##' @rdname k.qcs
##' @method k.qcs default
##' @inheritParams lab.qcdata
##' @param alpha The significance level (0.05 by default)
##' @param ... Other arguments passed to or from methods.
k.qcs.default <- function(x, var.index=1,replicate.index  =  2, material.index  =  3,
                          laboratory.index=4,  data.name = NULL, alpha = 0.05, ...)
{
  if (is.null(data.name)) data.name <- "Statistical Mandel k"

  obj<-lab.qcdata(data = x, var.index=var.index,replicate.index  =  replicate.index,
               material.index  =  material.index,
               laboratory.index=laboratory.index,  data.name = data.name)

  result<-k.qcs.lab.qcdata(x = obj,  alpha = alpha)

  return(result)
} #k.qcs

##' @rdname  k.qcs
##' @method k.qcs lab.qcdata
##' @inheritParams k.qcs.default
##' @export
k.qcs.lab.qcdata<- function(x, alpha = 0.05, ...)
  {

  if(is.null(x) || !inherits(x, "lab.qcdata"))
    stop("x must be an objects of class (or extending) 'lab.qcdata'")

  data.name <- attributes(x)$data.name
  x.lab.qcs <- lab.qcs(x)
  statistics <- x.lab.qcs$statistics.material
  s.i <- x.lab.qcs$statistics.Laboratory$s.i
  p <- x.lab.qcs$p
  n <- x.lab.qcs$n
  m <- x.lab.qcs$m

  v1<-(p-1)*(n-1)
  v2<-n-1

  kcrit <- sqrt(p/(1+(p-1)*qf(alpha,v1,v2,lower.tail=TRUE)))


  material <- row.names(x.lab.qcs$statistics.material)
  laboratory <- unique(x.lab.qcs[[1]]$laboratory)
    k.i<-matrix(,nrow =p ,ncol =m )
    for(i in 1:m)
    {
      ind <- x.lab.qcs$statistics.Laboratory$material==material[i]
      k.i[,i] <- s.i[ind]/statistics$S_r[i]
      }
    colnames(k.i) <- material
    row.names(k.i) <- laboratory
    violations <- k.i <= kcrit

    result <- list (lab.qcdata = x, lab.qcs = x.lab.qcs, p = p, n = n, m = m,
                    k = k.i, k.critical = kcrit, violations = violations, data.name = data.name )
    oldClass(result) <- c("lab.qcs")
    attr(result, "object.name") <- data.name
    attr(result, "type.data") <- "k.qcs"

    return(result)
  }

#-------------------------------------------------------------------------
# Cochran Test Statistic
#-------------------------------------------------------------------------
##' Function to compute the Cochran test statistic
##'
##' Function to estimate the Cochran test statistic.
##' @param x   An object of class \code{lab.qcdata} (Univariate Quality Control Data).
##' @export
##' @references
##' Wilrich Peter-T. (2013),  Critical values of mandel's h and k,
##' the grubbs and the Cochran test statistic. Asta-Advances in Statistical Analysis, 97(1):1-10.
##' 
##' ASTM E 691 (1999), Standard practice for conducting an interlaboratory study
##'   to determine the precision of a test method. American Society for Testing and Materials. West Conshohocken, PA, USA.
##' @examples
##'
##' library(ILS)
##' data(Glucose)
##' Glucose.qcdata <- lab.qcdata(Glucose)
##' str(Glucose.qcdata)
##' cochran.test(Glucose.qcdata)
cochran.test <- function(x, ...) {
  UseMethod("cochran.test")
}

##' @rdname cochran.test
##' @method cochran.test default
##' @inheritParams lab.qcdata
##' @param alpha The significance level (0.05 by default)
##' @param ... Other arguments passed to or from methods.
cochran.test.default <- function(x, var.index=1,replicate.index  =  2, material.index  =  3,
                          laboratory.index=4,  data.name = NULL, alpha = 0.05, ...)
{
  if (is.null(data.name)) data.name <- "Statistical Mandel k"

  obj<-lab.qcdata(data = x, var.index=var.index,replicate.index  =  replicate.index,
               material.index  =  material.index,
               laboratory.index=laboratory.index,  data.name = data.name)

  result<-cochran.test.lab.qcdata(x = obj,  alpha = alpha)

  return(result)
} #cochran.test

##' @rdname  cochran.test
##' @method cochran.test lab.qcdata
##' @inheritParams cochran.test.default
##' @export
cochran.test.lab.qcdata<-function(x, alpha = 0.05,...){

  if(!is.null(x) & !inherits(x, "lab.qcdata") & !is.list(x))
    stop("x must be an objects of class (or extending) 'lab.qcdata'")

  x.lab.qcs <- lab.qcs(x)
  stat <- x.lab.qcs$statistics.Laboratory
  material <- row.names(x.lab.qcs$statistics.material)
  laboratory <- unique(x$laboratory)

  S2max <- tapply(stat$s.i,stat$material,max)
  ind.max <- tapply(stat$s.i,stat$material,which.max)
  laboratory.max <- laboratory[ind.max]

  p <- x.lab.qcs$p
  n <- x.lab.qcs$n
  m <- x.lab.qcs$m

  C <- vector()
  p.value <- vector()

  v1 <- (p-1)*(n-1);
  v2 <- n-1
  Ccrit <- 1/(1+(p-1)*qf(alpha/p,v1,v2,lower.tail=TRUE))


  for(i in 1:m){
    C[i] <- S2max[i]/sum((stat$s.i[stat$material == material[i]])^2)
    p.value[i] <- round(pf(C[i],v1,v2,lower.tail=T),4)
  }

  result <- list(result = data.frame(Smax = laboratory.max, Material = material,
                       C = C,  p.value = p.value),C.critical = Ccrit, alpha.test = alpha/p)


  oldClass(result) <- c("cochran.test")

    return(result)
}

# @rdname lab.qcs
##' @method print cochran.test
## @param x A \code{test.cochran} object for which a print is desired.
##' @export

print.cochran.test <- function(x, ...) {
  cat("\nTest Cochran", "\n")
  cat("\n Critical value:",x[[2]],"\n")
  cat("\n Alpha test:",x[[3]],"\n")
  print(x[[1]])}

#-------------------------------------------------------------------------
# Grubbs Test Statistic
#-------------------------------------------------------------------------
##' Function to compute the Grubbs test statistic
##'
##' Function to estimate the Grubbs test statistic.
##' @param x   An object of class \code{lab.qcdata} (Univariate Quality Control Data).
##' @export
##' @references
##' Wilrich Peter-T. (2013), Critical values of Mandel's h and k,
##' the Grubbs and the Cochran test statistic. Asta-Advances in Statistical Analysis, 97(1):1-10.
##' 
##' ASTM E 691 (1999), Standard practice for conducting an interlaboratory study
##' to determine the precision of a test method. American Society for Testing and Materials. West Conshohocken, PA, USA.
##' @examples
##'
##' library(ILS)
##' data(Glucose)
##' Glucose.qcdata<- lab.qcdata(Glucose)
##' str(Glucose.qcdata)
##' grubbs.test(Glucose.qcdata)
grubbs.test <- function(x, ...) {
  UseMethod("grubbs.test")
}

##' @rdname grubbs.test
##' @method grubbs.test default
##' @inheritParams lab.qcdata
##' @param alpha The significance level (0.05 for default)
##' @param ... Other arguments passed to or from methods.
grubbs.test.default <- function(x, var.index=1,replicate.index  =  2, material.index  =  3,
                          laboratory.index=4,  data.name = NULL, alpha = 0.05, ...)
{
  if (is.null(data.name)) data.name <- "Statistical Mandel k"

  obj<-lab.qcdata(data = x, var.index=var.index,replicate.index  =  replicate.index,
               material.index  =  material.index,
               laboratory.index=laboratory.index,  data.name = data.name)

  result<-grubbs.test.lab.qcdata(x = obj,  alpha = alpha)

  return(result)
} #grubbs.test

##' @rdname  grubbs.test
##' @method grubbs.test lab.qcdata
##' @inheritParams grubbs.test.default
##' @export
grubbs.test.lab.qcdata <-function(x, alpha = 0.05,...){
  x.lab.qcs <- lab.qcs(x)
  stat <- x.lab.qcs$statistics.Laboratory
  material <- row.names(x.lab.qcs$statistics.material)
  laboratory <- unique(x$laboratory)

  p <- x.lab.qcs$p
  n <- x.lab.qcs$n
  m <- x.lab.qcs$m



  Gh <- vector()
  Gl <- vector()
  S <- vector()
  ph.value <- vector()
  pl.value <- vector()

  mean.i <- stat$mean.i
  mean <- x.lab.qcs$statistics.material$mean
  S <- x.lab.qcs$statistics.material$S



  ind.max <- tapply(stat$mean.i,stat$material,which.max)
  ind.min <- tapply(stat$mean.i,stat$material,which.min)
  laboratory.max <- laboratory[ind.max]
  laboratory.min <- laboratory[ind.min]


  for(i in 1:m){

    Gl[i] <- (mean[i] - mean.i[stat$material == material[i]][ind.min[i]])/S[i]
    pl.value[i] <- round(pt(Gl[i],(p-1),lower.tail=F),4)
    Gh[i] <- (mean.i[stat$material == material[i]][ind.max[i]] - mean[i] )/S[i]
    ph.value[i] <- round(pt(Gh[i],(p-1),lower.tail=F),4)
  }

  gcrit <- (n-1)*qt((1-alpha/p),(n-2))/sqrt(n*(n-2+(qt((1-alpha/p),(n-2)))^2))

  result <- list(result = data.frame(Material = material, Gmax = laboratory.max,
                                     G.max = Gh,
                                     p.value.max = ph.value, Gmin = laboratory.min,
                                     G.min = Gl,
                       p.value.min = pl.value),G.critical = gcrit,
                 alpha.test = alpha/p)

  oldClass(result) <- c("grubbs.test")
  return(result)
}


##' @method print grubbs.test
## @param x A \code{test.cochran} object for which a print is desired.
##' @export

print.grubbs.test <- function(x, ...) {
  cat("\nTest Grubbs", "\n")
  cat("\n Critical value:",x[[2]],"\n")
  cat("\n Alpha test:",x[[3]],"\n")
  print(x[[1]])}



#-------------------------------------------------------------------------
# AOV
#-------------------------------------------------------------------------
##' Function to compute the AOV
##'
##' Function to compute the analysis of variance of ILS data, taking into account the laboratories and material factors.
##' @param x An object of class \code{lab.qcdata} (Univariate Quality Control Data).
##' @export
##' @references
##' WHothorn T., Bretz, F., and Westfall, P. (2008), Simultaneous inference in general parametric models.
##' Biometrical Journal, 50(3):346-363.
##' 
##' Heyden, Y., Smeyers-Verbeke, J. (2007), Set-up and evaluation of interlaboratory studies. J. Chromatogr. A, 1158:158-167.
##' 
##' @examples
##' \dontrun{
##' library(ILS)
##' data(Glucose)
##' Glucose.qcdata <- lab.qcdata(Glucose)
##' str(Glucose.qcdata)
##' lab.aov(Glucose.qcdata,level = 0.95, plot = TRUE, pages = 1)
##' }

lab.aov <- function(x, ...) {
  UseMethod("lab.aov")
}
##' @rdname lab.aov
##' @method lab.aov default
##' @inheritParams lab.qcdata
##' @param level Requested confidence level (0.95 by default).
##' @param plot  If TRUE, confidence intervals are plot.
##' @param pages By default 0, it indicates the number of pages over which to spread the output. For example,
##' if pages=1,  all terms will be plotted on one page with the layout performed automatically.
##'  If pages=0, one plot will be displayed by each tested material.
##' @param ... Other arguments passed to or from methods.
lab.aov.default <- function(x, var.index=1,replicate.index  =  2, material.index  =  3,
                          laboratory.index=4,  data.name = NULL, level = 0.95,plot = FALSE, pages = 0, ...)
{
  if (is.null(data.name)) data.name <- "Statistical Mandel k"

  obj<-lab.qcdata(data = x, var.index=var.index,replicate.index  =  replicate.index,
               material.index  =  material.index,
               laboratory.index=laboratory.index,  data.name = data.name)

  result<-lab.aov.lab.qcdata(x = obj,  level = level,plot = plot, pages = pages)

  return(result)
} #lab.aov


##' @rdname lab.aov
##' @method lab.aov lab.qcdata
##' @inheritParams lab.aov.default
##' @export
lab.aov.lab.qcdata <- function(x,level = 0.95,plot = FALSE, pages = 0,...){

aovModel <- list()
conf <- list()
.Pairs <- list()
material <- unique(x$material)
m <- length(material)

if(plot ==TRUE){

  n.plots <- m
  if (pages > 0)
    if (pages > n.plots)
      pages <- n.plots
    if (pages < 0)
      pages <- 0
    if (pages != 0) {
      ppp <- n.plots%/%pages
      if (n.plots%%pages != 0) {
        ppp <- ppp + 1
        while (ppp * (pages - 1) >= n.plots) pages <- pages - 1
      }
      c <- r <- trunc(sqrt(ppp))
      if (c < 1)
        r <- c <- 1
      if (c * r < ppp)
        c <- c + 1
      if (c * r < ppp)
        r <- r + 1
      oldpar <- par(mfrow = c(r, c))
    }
    else {
      ppp <- 1
      oldpar <- par()
    }
}

for (i in 1:m){
  indm<-x$material==material[i]
  y <- x$x[indm]
  laboratory <- x$laboratory[indm]
  data <- data.frame(y,laboratory)

  aovModel[[i]] <- aov(y ~ laboratory,data=data)
  .Pairs[[i]] <- glht(aovModel[[i]], linfct = mcp(laboratory = "Tukey"))
  conf[[i]] <- confint(.Pairs[[i]],level = level) # confidence intervals
}


if(plot ==TRUE){
  old.oma <- par(oma=c(0,5,0,0))
  for (i in 1:m){
    title <- paste(level*100,"%"," ","Confidence Level",sep="")
    subtitle  = paste("Material",material[i])
    plot(confint(.Pairs[[i]],level = level), main=title,sub = subtitle)
  }
  par(old.oma)
}

par(mfrow=c(1,1))
  names(conf)<- paste("Material:",material)
  names(.Pairs)<-paste("Material:",material)
  names(aovModel)<-paste("Material:",material)
  for (i in 1:m) {cat("\n AOV of Material:",material[i])

    print(summary(aovModel[[i]]))
    print(summary(.Pairs[[i]])) # pairwise tests
    print(conf[[i]])
  }

    result <- list(Models = aovModel,Confidence =conf)

    return(result)
invisible()

}
