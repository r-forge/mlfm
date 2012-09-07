\name{mlfm}
\alias{mlfm}
\alias{print.mlfm}

\title{Semiparametric multilevel frailty models, estimated via EMPL method}

\description{The \code{mlfm} function fits proportional hazard models 
 with K>0 random effects.
 Inference is based on the EM-PPL algorithm proposed by Horny (2009).}

\usage{
  mlfm(formula, data, 
       eps = 1e-10, frailty.eps = 1e-11, maxit = 50,
       showtime = TRUE, verbose = FALSE, sparse = TRUE)
  
  \method{print}{mlfm}(x, digits = max(3, getOption("digits") - 3) - 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula object, with the response on the left of a ~ operator,                                      
    and the terms on the right. 
    The response must be a survival object as returned by the Surv() function
}
  \item{data}{a \code{\link{data.frame}} containing all the variables}
  \item{eps}{the threshold used to state convergence}
  \item{frailty.eps}{the threshold used to state convergence within PPL fits
    (see \code{\link[survival]{coxph.control}})}
  \item{maxit}{the maximum number of iteration for the outer loop (EM)}
  \item{showtime}{display the execution time}
  \item{verbose}{display the progression status over EM steps}
  \item{sparse}{wether to use sparse computation
    (see \code{\link[survival]{frailty}})}
  \item{x}{the object of class \code{mlfm} to print}
  \item{digits}{the number of digits to display in results}
  \item{...}{other arguments (see \code{\link{print}})}
}

\note{The standard errors of regression coefficients 
as given by the EMPL algorithm are underestimated, so are not returned.
A correct estimation is, for example, the one based on Louis (1982) approach,
still not implemented in \code{mlfm}.}

\references{Horny G (2009) 
\emph{Inference in mixed proportional hazard models with K random effects} 
Stat Papers \bold{50}:481--499}


\author{
Authors: Federico Rotolo and Guillaume Horny

Maintainer: Federico Rotolo <federico.rotolo@stat.unipd.it>
}

\keyword{EMPL}
\keyword{mixed effects}
\keyword{survival}
\keyword{frailty}
\keyword{multylevel}
\keyword{semiparametric}

\seealso{
  \code{\link[survival]{coxph}}, 
  \code{\link[parfm]{parfm}}, 
  \code{\link[frailtypack]{frailtypack}}}

\examples{
############################## - Data simulation - #############################
set.seed(1)
nsim <- 2000
data <- data.frame(id = 1:nsim)

### - Clustering - ###
grVar <- cbind(nclust = c(50, 20, 80),
               Fpar   = c(1, .5, 1.5))

for (gr in 1:nrow(grVar)) {
  data[[paste("gr", gr, sep="")]] <- apply(
    rmultinom(n= nsim, 1, rep(1, grVar[gr, 'nclust'])), 2, which.max)
  
  data[[paste("u", gr, sep="")]] <- rgamma(grVar[gr, 'nclust'],
                                           shape=1/grVar[gr, 'Fpar'],
                                           scale=grVar[gr, 'Fpar'])[
                                             data[[paste("gr", gr, sep="")]]]
}

head(data)

### - Covariates - ###
data$Treat <- rbinom(nsim, 1, .5)
data$Age   <- rnorm(nsim, 60, 4)

### - Times - ###
lambda <- 2
rho <- 1.2
tcens <- 1.5
simBetas <- c(.5, log(1.2) / 10)

data$Time <-pmin(tcens, (- log(runif(nsim, 0, 1)) / (lambda * 
  apply(data[, paste("u", 1:nrow(grVar), sep="")], 1, prod) *
  as.vector(exp(simBetas \%*\% t(data[, c("Treat", "Age")])))
                                                     ))^(1 / rho))
data$Status <- 0 + (data$Time < tcens)

summary(data)
head(data)

############################## - Model fitting - #############################
mlfmod <- mlfm(Surv(Time, Status) ~ Treat + Age +
  frailty(gr1) +frailty(gr2) +frailty(gr3) , data=data, verbose=TRUE)
mlfmod
}
