################################################################################
#  Print method for class 'mlfm'                                               #
################################################################################
#                                                                              #
#  This function prints the objects of class 'mlfm'                            #
#                                                                              #
#  Its parameters are                                                          #
#   - x         : the fitted model, object of class 'parfm'                    #
#   - digits    : number of significant digits                                 #
#                                                                              #
################################################################################
#   Author: Federico Rotolo <federico.rotolo@stat.unipd.it>                    #
#                                                                              #
#   Date: April 17, 2012                                                       #
#   Last modification on: April 27, 2012                                       #
################################################################################

print.mlfm <- function(x, 
                       digits = max(3, getOption("digits") - 3),
                       ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  print.default(cbind(
    coef = format(x$beta, digits = digits),
    se   = NA,
    p    = NA
    ), na.print = "", print.gap = 2,  quote = FALSE)
  
  cat(paste("\nIterations:",
            paste(x$iter, names(x$iter), collapse=", ")))
  if (length(x$theta) - 1) {
    cat("\n     Variances of random effects")
    names(x$theta) <- paste("     ", names(x$theta))
    print.default(cbind(" "=
      format(x$theta, digits = digits)
                        ), na.print="", print.gap = 2,  quote = FALSE)
  } else {
    cat(paste("\n     Variance of random effect=", 
              format(x$theta, digits = digits)))
  }
}
