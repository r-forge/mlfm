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
#   Last modification on: August 28, 2012                                      #
################################################################################

print.mlfm <- function(x, 
                       digits = max(3, getOption("digits") - 3),
                       ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  print.default(cbind(
    coef = signif(x$beta, digits = digits),
    se   = signif(x$ses,  digits = digits),
    p    = signif(pchisq((x$beta / x$ses)^2, df=1, lower.tail=FALSE), 
                  digits = digits)
    ), na.print = "", print.gap = 2,  quote = FALSE)
  
  cat(paste("\nIterations:",
            paste(x$iter, names(x$iter), collapse=", ")))
  if (length(x$theta) - 1) {
    cat("\n     Variances of random effects")
    names(x$theta) <- paste("     ", names(x$theta))
    print.default(cbind(" "=
      format(x$theta, digits = digits)
                        ), na.print="", print.gap = 2,  quote = FALSE)
    cat(paste("Total frailty variance:",
              format(prod(x$theta + 1) - 1, digits = digits), "\n"))
  } else {
    cat(paste("\n     Variance of random effect=", 
              format(x$theta, digits = digits)))
  }
  
  if (x$convergence)
    cat("\nNOTE: the estimation procedure did not converge!\n")
}
