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
##' @param x  Object lab.qcd (Functional Quality Control Data)
##' @param ... Arguments passed to or from methods.
##' @export
##' @examples
##' 
##' library(ILS)
##' data(Glucose)
##' Glucose.qcd <- lab.qcd(Glucose)
##' str(Glucose.qcd)
##' Glucose.qcs <- lab.qcd(Glucose.qcd)
##' str(Glucose.qcs)
##' summary(Glucose.qcs)
lab.qcs <- function(x, ...)
  #.........................................................................  
  {
  if(is.null(x) || !inherits(x, "lab.qcd"))
    stop("x must be an objects of class (or extending) 'lab.qcd'")
  
  p <- length(unique(x$laboratory))
  m <- length(unique(x$material))
  n <- length(unique(x$replicate))
  material<-unique(x$material)
  laboratory<-unique(x$laboratory)
  mean.i<-matrix(,nrow =p ,ncol =m )
  s.i<-matrix(,nrow = p,ncol = m)
  S_r<-vector()
  S_R<-vector()
 
   for(i in 1:m)
  {
    ind<-x$material==material[i]
   mean.i[,i]<- tapply(x$x[ind],x$laboratory[ind],mean)
  s.i[,i] <- tapply(x$x[ind],x$laboratory[ind],sd)
  S_r[i]<-sqrt(mean(s.i[,i]^2))
  S_R[i]<-sqrt(var(mean.i[,i])+((n-1)/n)*S_r[i]^2)
  
  }
 
 colnames(mean.i) <- material
 row.names(mean.i) <- laboratory
 colnames(s.i) <- material
 rownames(s.i) <-laboratory
 statistics <- data.frame(mean=apply(mean.i,2,mean),S=apply(mean.i,2,sd),
                       S_r=S_r, S_R=S_R)
 statisticsL <- data.frame(mean = mean.i,sd = s.i)
 result <- list (lab.qcd = x, statistics =  statistics, p = p, n = n, m = m,
                 mean.i = mean.i, S.i = s.i, statistics.Laboratory = statisticsL  ) 
 oldClass(result)<-c("lab.qcs")
 attr(result, "object.name") <- "lab.qcs"
 attr(result, "type.data") <- "univariate"
 
return(result)
} # lab.qcs
#.........................................................................

##' @rdname lab.qcs
##' @method print lab.qcs
## @param x A \code{lab.qcd} object for which a print is desired.
##' @export
print.lab.qcs <- function(x, ...) str(x,1)
#.........................................................................
##' @rdname lab.qcs
##' @method summary lab.qcs
##' @param object A \code{lab.qcs} object for which a summary is desired.
##' @export
summary.lab.qcs <- function(object, ...)
  #.........................................................................
{
 
  
  type.data <- attributes(object)$object.name

  cat("\nSummary of material statistics:\n")
  print(object$statistics)
  cat("\nNumber of laboratories: ", object[[3]])
  cat("\nNumber of materials: ", object[[5]])
  cat("\nNumber of replicate: ", object[[4]])
  
  result <- switch(type.data, 
                  "lab.qcs" =  {
                    cat("\nSummary for Laboratory:\n")
                    print(object[[8]])
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
##' @param x   R object (used to select the method). See details.
##' @export
##' @references 
##' \describe{
##'   \item{}{Wilrich Peter-T. (2013),  Critical values of Mandel's h and k, 
##'   the Grubbs and the Cochram test statistic. Asta-Advances in Statistical Analysis, 97(1):1-10.}
##'   \item{}{ASTM E 691 (1999), Standard practice for conducting an interlaboratory study 
##'   to determine the precision of a test method. American Society for Testing and Materials. West Conshohocken, PA, USA.}
##' }
##' @examples
##' 
##' library(ILS)
##' data(Glucose)
##' Glucose.qcd <- lab.qcd(Glucose)
##' str(Glucose.qcd)
##' h<- h.qcs(Glucose.qcd, alpha = 0.005)
##' summary(h)
##' plot(h)

h.qcs <- function(x, ...) {
  UseMethod("h.qcs")
}

##' @rdname h.qcs
##' @method h.qcs default
##' @inheritParams lab.qcd
##' @param alpha The significance level (0.05 by default)
##' @param ... Arguments passed to or from methods.
h.qcs.default <- function(x, var.index=1,replicate.index  =  2, material.index  =  3,
                          laboratory.index=4,  data.name = NULL, alpha = 0.05, ...)
  {
  if (is.null(data.name)) data.name <- "Statistical Mandel h"

    obj<-lab.qcd(data = x, var.index=var.index,replicate.index  =  replicate.index, 
                 material.index  =  material.index,
                 laboratory.index=laboratory.index,  data.name = data.name)
    
    result<-h.qcs.lab.qcd(x = obj,  alpha = alpha)
    
  return(result)
} #h.qcs

##' @rdname  h.qcs
##' @method h.qcs lab.qcd
##' @inheritParams h.qcs.default
##' @export
h.qcs.lab.qcd <- function(x, alpha = 0.05, ...)
{
  data.name <- attributes(x)$data.name
  x.lab.qcs <- lab.qcs(x)
  statistics <- x.lab.qcs[[2]]
  mean.i <- x.lab.qcs[[6]]
  p <- x.lab.qcs[[3]]
  n <- x.lab.qcs[[4]]
  m <- x.lab.qcs[[5]]
  hcrit <- (p-1)*qt((1-alpha/2),(p-2))/sqrt(p*(p-2+qt((1-alpha/2),(p-2))^2))
  material <- row.names(x.lab.qcs[[2]])
  laboratory <- unique(x.lab.qcs[[1]]$laboratory)
  h.i <- matrix(,nrow = p,ncol = m)
  for(i in 1:m)
  {
    h.i[,i] <- (mean.i[,i]-statistics[i,1])/statistics[i,2]
  }
  
  i <- 1
  (mean.i[,i]-statistics[i,1])/statistics[i,2]
  
  
  colnames(h.i) <- material
  rownames(h.i) <- laboratory
  violations <- h.i<=hcrit
  result <- list (lab.qcd = x, lab.qcs = x.lab.qcs, p = p, n = n, m = m,
                  h = h.i, h.critial = hcrit, violations = violations, data.name = data.name ) 
  
  oldClass(result) <- c("lab.qcs")
  attr(result, "object.name") <- "h.qcs"
  attr(result, "type.data") <- "univariate"
  
  
  return(result)
  
} #h.qcs

# Statistic k
#-------------------------------------------------------------------------
##' Function to calcute the Mandel's k statistic
##'
##' This function is used to compute the statistic k of Mandel.
##' @param x   an R object (used to select the method). See details.
##' @export
##' @references 
##' \describe{
##'   \item{}{Wilrich Peter-T. (2013),  Critical values of Mandel's h and k, 
##'   the Grubbs and the Cochram test statistic. Asta-Advances in Statistical Analysis, 97(1):1-10.}
##'   \item{}{ASTM E 691 (1999), Standard practice for conducting an interlaboratory study 
##'   to determine the precision of a test method. American Society for Testing and Materials. West Conshohocken, PA, USA.}
##' }
##' @examples
##' 
##' library(ILS)
##' data(Glucose)
##' Glucose.qcd <- lab.qcd(Glucose)
##' str(Glucose.qcd)
##' k<- k.qcs(Glucose.qcd, alpha = 0.005)
##' summary(k)
##' plot(k)
k.qcs <- function(x, ...) {
  UseMethod("k.qcs")
}

##' @rdname k.qcs
##' @method k.qcs default
##' @inheritParams lab.qcd
##' @param alpha The significance level (0.05 by default)
##' @param ... arguments passed to or from methods.
k.qcs.default <- function(x, var.index=1,replicate.index  =  2, material.index  =  3,
                          laboratory.index=4,  data.name = NULL, alpha = 0.05, ...)
{
  if (is.null(data.name)) data.name <- "Statistical Mandel k"
  
  obj<-lab.qcd(data = x, var.index=var.index,replicate.index  =  replicate.index, 
               material.index  =  material.index,
               laboratory.index=laboratory.index,  data.name = data.name)
  
  result<-k.qcs.lab.qcd(x = obj,  alpha = alpha)
  
  return(result)
} #k.qcs

##' @rdname  k.qcs
##' @method k.qcs lab.qcd
##' @inheritParams k.qcs.default
##' @export
k.qcs.lab.qcd<- function(x, alpha = 0.05, ...)
  {
  data.name <- attributes(x)$data.name
  x.lab.qcs <- lab.qcs(x)
  statistics <- x.lab.qcs[[2]]
  s.i <- x.lab.qcs[[7]]
  p <- x.lab.qcs[[3]]
  n <- x.lab.qcs[[4]]
  m <- x.lab.qcs[[5]]
  
  v1<-(p-1)*(n-1)
  v2<-n-1
  
  kcrit<-sqrt(p/(1+(p-1)*qf(alpha,v1,v2,lower.tail=TRUE)))

  
  material <- row.names(x.lab.qcs[[2]])
  laboratory <- unique(x.lab.qcs[[1]]$laboratory)
    k.i<-matrix(,nrow =p ,ncol =m )
       for(i in 1:m)
    {
      ind<-x$material==material[i]
      k.i[,i]<-s.i[,i]/statistics$S_r[i]
      }
    colnames(k.i)<-material
    row.names(k.i)<-laboratory
    violations <- k.i<=kcrit
    
    result <- list (lab.qcd = x, lab.qcs = x.lab.qcs, p = p, n = n, m = m,
                    k = k.i, k.critical = kcrit, violations = violations, data.name = data.name ) 
    oldClass(result) <- c("lab.qcs")
    attr(result, "object.name") <- "k.qcs"
    attr(result, "type.data") <- "univariate"
    
    return(result)
  }

#-------------------------------------------------------------------------
# Cochram Test Statistic
#-------------------------------------------------------------------------
##' Function to compute the Grubbs test statistic.
##'
##' Function to estimate the Cochram test statistic.
##' @param x   R object (used to select the method). See details.
##' @export
##' @references 
##' \describe{
##'   \item{}{Wilrich Peter-T. (2013),  Critical values of mandel's h and k, 
##'   the grubbs and the cochram test statistic. Asta-Advances in Statistical Analysis, 97(1):1-10.}
##'   \item{}{ASTM E 691 (1999), Standard practice for conducting an interlaboratory study 
##'   to determine the precision of a test method. American Society for Testing and Materials. West Conshohocken, PA, USA.}
##' } 
##' @examples
##' 
##' library(ILS)
##' data(Glucose)
##' Glucose.qcd <- lab.qcd(Glucose)
##' str(Glucose.qcd)
##' Cochram.test(Glucose.qcd)
Cochram.test <- function(x, ...) {
  UseMethod("Cochram.test")
}

##' @rdname Cochram.test
##' @method Cochram.test default
##' @inheritParams lab.qcd
##' @param alpha The significance level (0.05 by default)
##' @param ... Arguments passed to or from methods.
Cochram.test.default <- function(x, var.index=1,replicate.index  =  2, material.index  =  3,
                          laboratory.index=4,  data.name = NULL, alpha = 0.05, ...)
{
  if (is.null(data.name)) data.name <- "Statistical Mandel k"
  
  obj<-lab.qcd(data = x, var.index=var.index,replicate.index  =  replicate.index, 
               material.index  =  material.index,
               laboratory.index=laboratory.index,  data.name = data.name)
  
  result<-Cochram.test.lab.qcd(x = obj,  alpha = alpha)
  
  return(result)
} #Cochram.test

##' @rdname  Cochram.test
##' @method Cochram.test lab.qcd
##' @inheritParams Cochram.test.default
##' @export
Cochram.test.lab.qcd<-function(x, alpha = 0.05,...){

  x.lab.qcs <- lab.qcs(x)
  material<-row.names(x.lab.qcs$statistics)
  
  S.i2 <- x.lab.qcs[[7]]^2
  p <- x.lab.qcs[[3]]
  n <- x.lab.qcs[[4]]
  m <- x.lab.qcs[[5]]
  
  S2max <- vector()
  C <- vector()
  p.value <- vector()

  v1 <- (p-1)*(n-1);
  v2 <- n-1
  Ccrit <- 1/(1+(p-1)*qf(alpha,v1,v2,lower.tail=TRUE))  
  
  for(i in 1:m){
    S2max[i] <- max(S.i2[,i])
    C[i] <- S2max[i]/sum(S.i2[,i])
    p.value[i] <- round(pf(C[i],v1,v2,lower.tail=T),4)
  }

  result <- data.frame(Material = material,C = C, C.critical = Ccrit, p.value = p.value)

  cat("\nTest Cochram:", "\n")
    return(result)
}
#-------------------------------------------------------------------------
# Grubbs Test Statistic
#-------------------------------------------------------------------------
##' Function to compute the Grubbs test statistic.
##' 
##' Function to estimate the Grubbs test statistic.
##' @param x   an R object (used to select the method). See details.
##' @export
##' @references 
##' \describe{
##'   \item{}{Wilrich Peter-T. (2013), Critical values of Mandel's h and k, 
##'   the Grubbs and the Cochram test statistic. Asta-Advances in Statistical Analysis, 97(1):1-10.}
##'   \item{}{ASTM E 691 (1999), Standard practice for conducting an interlaboratory study 
##'   to determine the precision of a test method. American Society for Testing and Materials. West Conshohocken, PA, USA.}
##' }
##' @examples
##' 
##' library(ILS)
##' data(Glucose)
##' Glucose.qcd <- lab.qcd(Glucose)
##' str(Glucose.qcd)
##' Grubbs.test(Glucose.qcd)
Grubbs.test <- function(x, ...) {
  UseMethod("Grubbs.test")
}

##' @rdname Grubbs.test
##' @method Grubbs.test default
##' @inheritParams lab.qcd
##' @param alpha The significance level (0.05 for default)
##' @param ... arguments passed to or from methods.
Grubbs.test.default <- function(x, var.index=1,replicate.index  =  2, material.index  =  3,
                          laboratory.index=4,  data.name = NULL, alpha = 0.05, ...)
{
  if (is.null(data.name)) data.name <- "Statistical Mandel k"
  
  obj<-lab.qcd(data = x, var.index=var.index,replicate.index  =  replicate.index, 
               material.index  =  material.index,
               laboratory.index=laboratory.index,  data.name = data.name)
  
  result<-Grubbs.test.lab.qcd(x = obj,  alpha = alpha)
  
  return(result)
} #Grubbs.test
##' @rdname  Grubbs.test
##' @method Grubbs.test lab.qcd
##' @inheritParams Grubbs.test.default
##' @export
Grubbs.test.lab.qcd <-function(x, alpha = 0.05,...){
  x.lab.qcs <- lab.qcs(x)
  material <- row.names(x.lab.qcs$statistics)

  p <- x.lab.qcs[[3]]
  n <- x.lab.qcs[[4]]
  m <- x.lab.qcs[[5]]

  mean_max <- vector()
  mean_min <- vector()
  mean_mean <- vector()
  
  Gh <- vector()
  Gl <- vector()
  S <- vector()
  ph.value <- vector()
  pl.value <- vector()
  
  for(i in 1:m){
    mean_mean[i] <- mean(x.lab.qcs$mean.i[,i])
    S[i] <- sd(x.lab.qcs$mean.i[,i])
    mean_min[i] <- min(x.lab.qcs$mean.i[,i])
    Gl[i] <- (mean_mean[i] - mean_min[i])/S[i]
    pl.value[i] <- round(pt(Gl[i],(p-1),lower.tail=F),4)
    mean_max[i] <- max(x.lab.qcs$mean.i[,i])
    Gh[i] <- (mean_max[i]-mean_mean[i])/S[i]
    ph.value[i] <- round(pt(Gh[i],(p-1),lower.tail=F),4)
  }
  
  gcrit <- (n-1)*qt((1-alpha/p),(n-2))/sqrt(n*(n-2+(qt((1-alpha/p),(n-2)))^2))

  result <- data.frame(Material = material, G.max = Gh,p.value.max = ph.value,G.min = Gl,
                       p.value.min = pl.value,G.critical = gcrit)
  cat("\nTest Grubbs:", "\n")
  return(result)
}

#-------------------------------------------------------------------------
# AOV
#-------------------------------------------------------------------------
##' Function to compute the AOV
##' 
##' Function to compute the analysis of variance of ILS data, taking into account the laboratories and material factors.
##' @param x Object lab.qcd.
##' @export
##' @references 
##' \describe{
##'   \item{}{WHothorn T., Bretz, F., and Westfall, P. (2008), Simultaneous inference in general parametric models. 
##'   Biometrical Journal, 50(3):346-363.}
##'   \item{}{Heyden, Y., Smeyers-Verbeke, J. (2007), Set-up and evaluation of interlaboratory studies. J. Chromatogr. A, 1158:158-167.}
##' }
##' @examples
##' \dontrun{
##' library(ILS)
##' data(Glucose)
##' Glucose.qcd <- lab.qcd(Glucose)
##' str(Glucose.qcd)
##' lab.aov(Glucose.qcd,level = 0.95, plot = TRUE, pages = 1)
##' }
lab.aov <- function(x, ...) {
  UseMethod("lab.aov")
}
##' @rdname lab.aov
##' @method lab.aov default
##' @inheritParams lab.qcd
##' @param level Requested confidence level (0.95 by default)
##' @param plot  If TRUE, confidence intervals are plot.
##' @param pages By default 0, it indicates the number of pages over which to spread the output. For example, 
##' if pages=1,  all terms will be plotted on one page with the layout performed automatically.
##'  If pages=0, one plot will be displayed by each tested material.
##' @param ... Arguments passed to or from methods.
lab.aov.default <- function(x, var.index=1,replicate.index  =  2, material.index  =  3,
                          laboratory.index=4,  data.name = NULL, level = 0.95,plot = FALSE, pages = 0, ...)
{
  if (is.null(data.name)) data.name <- "Statistical Mandel k"
  
  obj<-lab.qcd(data = x, var.index=var.index,replicate.index  =  replicate.index, 
               material.index  =  material.index,
               laboratory.index=laboratory.index,  data.name = data.name)
  
  result<-lab.aov.lab.qcd(x = obj,  level = level,plot = plot, pages = pages)
  
  return(result)
} #lab.aov


##' @rdname lab.aov
##' @method lab.aov lab.qcd
##' @inheritParams lab.aov.default
##' @export
lab.aov.lab.qcd <- function(x,level = 0.95,plot = FALSE, pages = 0,...){

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
