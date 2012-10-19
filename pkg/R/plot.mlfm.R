################################################################################
#  Plot results of semiparametric multi-level frailty models                   #
################################################################################
#                                                                              #
#  This function plots results of multi-level frailty models                   #
#                                                                              #
#  Its parameters are                                                          #
#   -                                                                          #
#                                                                              #
################################################################################
#   Author: Federico Rotolo <federico.rotolo@stat.unipd.it>                    #
#   Original code by Guillaume Horny                                           #
#                                                                              #
#   Date: October 19, 2012                                                     #
#   Last modification on: October 19, 2012                                     #
################################################################################

plot.mlfm <- function(model,
                       xlim=NULL, 
                       main=NULL, 
                       coef=NULL, 
                       signcol=TRUE, 
                       cex=1){
  if (is.null(coef)) {
    n <- length(model$beta)
    if ("survreg"%in%class(model)) n<-n-1
    coef <- 1:n
  } else n <-length(coef)
  coef <- coef[n:1]
  intervals <- intervals(model)
  intervals <- intervals[coef,-2]
  if (is.null(xlim)){
    range <- c(0, max(max(intervals[is.finite(intervals)], na.rm=TRUE)*1.2, 2))
    if (!is.finite(range[2])) range <- c(0,100)
  }else{
    range <-xlim
  }
  if (is.null(main)) main<-""
  
  if (signcol) {color <- apply(intervals[,c("lower .95", "upper .95")], 1, function(x){
    if (is.na(x[1]) | is.na(x[2]) | (x[1]<1 & x[2]>1))  "grey" else "black"
  })} else color <- "black"
  
  par(mar=c(5,5,4,2)+.1)
  plot(as.vector(intervals), rep(1:n, 3),
       xlim=range, 
       ylim=.5+c(0, n),
       bty="]",
       xlab="HR", 
       ylab="", 
       yaxt="n",
       main=main,
       pch= c(rep(20, n), rep(91, n), rep(93, n)),
       cex=cex,
       col=color
       )
  abline(h=1:n, v=0:1, lty=c(rep(3,n),1,2), col=c(rep("grey",n),rep("black",2)))
  intervals[!is.finite(intervals[,3]),3] <- 
    intervals[!is.finite(intervals[,3]),1]+10^10
  segments(intervals[,2], 1:n, intervals[,3], 1:n, col=color)
  mtext(substr(dimnames(intervals)[[1]], 1, 10),
        side=2,
        at=1:n,
        las=1,
        cex=.7,
        col=color
        )	
}

intervals <- function(model){
  if ("mlfm" %in% class(model)) {
    intervals <- exp(model$beta + 
      model$ses %*% 
      t(qnorm(c(.5, .5, "lower .95"=.025, "upper .95"=.975)) ))
    row.names(intervals) <- names(model$beta)
    return(intervals)
  }
}
