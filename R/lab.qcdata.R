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

#
#  Main function to create a 'lab.qcdata' object
#
##' Quality Control Data
##' 
##' It creates a 'lab.qcdata' class object to perform the interlaboratory study. 
##' This object is used to plot ILS data and more.
##' 
##' @aliases lab.qcdata 
##' @param data Matrix or data-frame that contains the data, replicate index, type of material, and the laboratory.
##' @param var.index Scalar with the column number corresponding to the observed variable (the critical to quality variable). 
##' Alternatively, a string with the name of a quality variable can be provided.
##' @param replicate.index Scalar with the column number corresponding to the index each replicate.
##' @param material.index Scalar corresponding to the replicated number.
##' @param laboratory.index Scalar that defines the index number of each laboratory.
##' @param data.name String specifying the name of the variable which appears on the plots. 
##' If name is not provided, it is retrieved from the object.
##' @export
##' @examples
##' library(ILS)
##' data(Glucose)
##' Glucose.qcdata <- lab.qcdata(Glucose)
##' str(Glucose.qcdata)
##' summary(Glucose.qcdata)

lab.qcdata <- function(data, var.index=1,replicate.index  =  2, material.index  =  3,
                    laboratory.index=4,  data.name = NULL) 
  #.........................................................................
{
  
  if (!is.matrix(data) & !is.data.frame(data))
    stop("object must be a matrix or data.frame")
  
     result <- data[, c(var.index,replicate.index,material.index,
                        laboratory.index)]   
     
     names(result) <- c("x", "replicate","material","laboratory") 
     
     result$replicate <- as.factor(result$replicate)
     result$material <- as.factor(result$material)
     result$laboratory <- as.factor(result$laboratory)
     

  
  if (is.null(data.name))
    data.name <- deparse(substitute(data))
  
  attr(result, "data.name") <- data.name
  
  oldClass(result) <- c("lab.qcdata", "data.frame") #cambie la clase del resultado.
  
  return(result)
} # lab.qcdata
#.........................................................................
