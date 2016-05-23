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
# plot.lab.qcs
#-------------------------------------------------------------------------
##' Plot method for 'lab.qcs' objects
##' 
##' Generic function for plotting objects of 'lab.qcs' class. Results of  univariate ILS studies are graphically shown.
##' 
##' @method plot lab.qcs
##' @param x  Functional data object or a list with functional objects 
##' @param title Main title for the plot
##' @param xlab Title for the x axis
##' @param ylab Title for the y axis
##' @param col Color specifications
##' @param ylim Numeric vectors of length 2 (coordinates ranges)
##' @param ...  arguments to be passed to or from methods.
##' @export
plot.lab.qcs <- function(x, title = NULL, xlab = NULL, ylab = NULL, col = NULL, ylim = NULL, ...)
  #..............................................................................      
{
  if(!is.null(x) & !inherits(x, "h.qcs") & !inherits(x, "k.qcs") & !is.list(x))
    stop("x must be an objects of class (or extending) 'qcs'")

oldpar <- par(mar = c(4, 3, 1, 1) + 0.1)
data.name <- x[[9]]
object.name <- attributes(x)$object.name 
if (is.null(title)) title <- data.name
if (is.null(xlab)) xlab <- "Laboratory"
if (is.null(ylab)) ylab <- "Statistical"


st<-t(x[[6]])
p <-x$p
m <- x$m
crit <- x[[7]]


if (is.null(col)){
  if (m==1){
    col <- terrain.colors(p)
  }else{
    col <- terrain.colors(m)
  }
}

legend.text = rownames(t(x[[6]]))

if (object.name=="h.qcs")
{
  if (is.null(ylim)) ylim <- c(-3,3)
if (m>1){
  
  
  barplot(height = st,width = 0.05,main=title,names.arg = colnames(t(x[[6]])),
        beside = T,
        col=col,xlab=xlab,
        ylim=ylim,axisnames=T)

  legend ("topleft",legend = legend.text,bty = "n",pch = 22,pt.bg = col,cex = 0.8)
}else{
  
  barplot(height = st,width = 0.05,main=title,names.arg = colnames(t(x[[6]])),
          beside = T,
          col=col,xlab=xlab,
          ylim=ylim,axisnames=T)
  
}

abline(h=c(crit,-crit),lty="dashed")

}else{

  if (is.null(ylim)) ylim <- c(0,3)
  if (m>1){
  barplot(height = st,width = 0.05,main=title,names.arg = colnames(t(x[[6]])),
        beside = T,
        horiz=F,col=col,xlab=xlab,
        ylim=ylim,axisnames=T)
    legend ("topleft",legend = legend.text,bty = "n",pch = 22,pt.bg = col,cex = 0.8,pt.cex = 1.5)
  }else{
    barplot(height = st,width = 0.05,main=title,names.arg = colnames(t(x[[6]])),
            beside = T,
            horiz=F,col=col,xlab=xlab,
            ylim=ylim,axisnames=T)
    
  }  

abline(h=crit,lty="dashed")
}



par(oldpar)

  
  #.........................................................................
} # plot.lab.qcs
#...........................................................................
