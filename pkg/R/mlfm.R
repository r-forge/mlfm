################################################################################
#  Semiparametric multi-level frailty models                                   #
################################################################################
#                                                                              #
#  This function implements the semiparametric EMPL estimation method          #
#  for multi-level frailty models                                              #
#                                                                              #
#  Its parameters are                                                          #
#   - formula  : a formula object, with the response                           #
#                on the left of a ~ operator,                                  #
#                and the terms on the right.                                   #
#                The response must be a survival object                        #
#                as returned by the Surv() function;                           #
#   - data     : a data.frame containing all the variables;                    #
#   - maxit    : the maximum number of iteration for the EM procedure          #
#   - showtime : display the execution time                                    #
#   - verbose  : display the progression status over EM steps                  #
#                                                                              #
################################################################################
#   Author: Federico Rotolo <federico.rotolo@stat.unipd.it>                    #
#   Original code by Guillaume Horny                                           #
#                                                                              #
#   Date: April 17, 2012                                                       #
#   Last modification on: July 9, 2012                                         #
################################################################################

mlfm <- function(formula, 
                 data, 
                 eps      = 1E-10,
                 maxit    = 50,
                 showtime = TRUE,
                 verbose  = FALSE) {
  require(survival)
  
  exTimeStart <- Sys.time()
  
  if (missing(data)) {
    data <- eval(parse(text=paste("data.frame(", 
                                  paste(all.vars(formula), collapse=", "),
                                  ")")))
  }
  
  ### - Checks - ###############################################################
  Call <- match.call()
  if (!match("formula", names(Call), nomatch=0))
    stop("A formula argument is required")
  
  Terms <- terms(formula, specials=c("frailty"), data = data)  
  
  families <- attr(Terms, "specials")$frailty
  if (is.null(families))
    stop("At least one group variable is required")
  ######################################################## - End of Checks - ###
  
  
  ### - Initialisation - #######################################################
  beta <- coefficients(coxph(Terms[-(families-1)], data = data))
  
  lnv <- mapply(rep, 0, apply(
    data[, substr(untangle.specials(Terms, "frailty")[[1]], 9, 
                  nchar(untangle.specials(Terms, "frailty")[[1]])-1)],
    2, function(x) nlevels(as.factor(x))))
  names(lnv) <-
    substr(untangle.specials(Terms, "frailty")[[1]], 9, 
           nchar(untangle.specials(Terms, "frailty")[[1]])-1)
  
  theta <- rep(NA, length(families))
  names(theta) <- names(lnv)
  for (f in names(lnv)) {
    sfmod <- coxph(update.formula(Terms[-(families-1)], paste(
      ".~. + frailty.gamma(", f, ", sparse=TRUE)")), data = data)
    
    theta[which(names(lnv) == f)] <- sfmod$history[[1]][[3]][sfmod$iter[1], 1]
  }  
  
  betaOld <- beta
  thetaOld <- theta
  lnvOld <- lnv
  
  iter <- c(EM = 0, outer = 0, 'Newton-Raphson' = 0)
  crit <- eps + .1
  ################################################ - End of Initialisation - ###
  
  if (verbose) cat("\n# - CRITERION - #\n")
  
  ### - EM cycle - #############################################################
  # corresponding to Step D of Horny (2009)
  while ((crit >= eps) && (iter[1] < maxit)) {
    gc()
    iter[1] <- iter[1] + 1 
    if (verbose)
      cat(paste(" Iteration ", format(iter[1], width=3), ": ", sep=""))
    # backup of the present betas
    betaOld <- beta
    
    ### - EXP step for v, each level separately - ##############################
    # corresponding to Steps A and B of Horny (2009)
    for (f in names(lnv)) {
      gc()
      # backup of the present log-frailties
      lnvOld[[f]] <- lnv[[f]]
      # backup of the present thetas
      thetaOld[f] <- theta[f]
      
      # PPL for present family
      data$offs <- 0
      for (of in (names(lnv)[names(lnv) != f])) {
        data$offs[as.numeric(data[, of]) > 0] <- 
          data$offs[as.numeric(data[, of]) > 0] + lnv[[of]][data[, of]]
      }
      
      form <- update.formula(formula, paste(
        ".~. + frailty.gamma(", f, ", sparse=TRUE, eps=1e-20)", 
        paste(" - frailty(", names(lnv), ")", 
              sep="", collapse=""),
        " + offset(offs)", sep="", collapse=""))
      
      mod <- coxph(form, data=data)
      
      iter[-1] <- iter[-1] + mod$iter
      
      # set-up of the new log-frailties and theta
      lnv[[f]] <- mod$frail
      theta[f] <- mod$history[[1]][[3]][mod$iter[1], 1]
    }
    #################################################### - End of EXP step - ###
    
    ### - MAX step for betas - #################################################
    # corresponding to Step C of Horny (2009)
    data$offs <- 0
    for (of in names(lnv)) {
      data$offs <- data$offs + lnv[[of]][data[, of]]
    }
    
    form <- update.formula(formula, paste(
      ".~.",
      paste(" - frailty(", names(lnv), ")", 
            sep="", collapse=""),
      " + offset(offs)", sep="", collapse=""))
    
    mod <- coxph(form, data=data, iter.max = 30, outer.max=20)
    
    iter[-1] <- iter[-1] + mod$iter
    
    # set-up of the new betas
    beta <- coefficients(mod)
    #################################################### - End of MAX step - ###
    
    # Termination criteria
    crit <- sum((c(log(theta), beta, unlist(lnv)) - 
      c(log(thetaOld), betaOld, unlist(lnvOld)))^2)
    if (verbose) cat(formatC(crit, format="E", digits=2), "\n", sep="")
  }
  ###################################################### - End of EM cycle - ###
  if (verbose && (iter[1] >= maxit))
    cat("\nThe procdeure did not converge!\n")
  
  mod <- list(
    call        = Call,
    beta        = coefficients(mod),
    theta       = theta,
    frail       = lnv,
    iter        = iter,
    convergence = as.numeric(iter[1] >= maxit),
    crit        = crit,
    extime      = diff(c(exTimeStart, Sys.time()))
    )
  class(mod) <- "mlfm"
  
  if (showtime) cat(paste("\nExecution time:", 
                          format(mod$extime, 
                                 digits = max(3, getOption("digits") - 3)), 
                          "\n"))
  
  ### - Return results - #######################################################
  return(mod)
}
