\name{ARCensReg}
\alias{ARCensReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Censored Linear Regression Model with Autoregressive Errors
}
\description{
It fits an univariate left or right censored linear regression model with autoregressive errors under the normal distribution using the SAEM algorithm. It provides estimates and standard errors of the parameters, prediction of future observations and it supports missing values on the dependent variable. It also provides convergence plots when exists at least one censored observation.
}
\usage{
ARCensReg(cc, y, x, p = 1, cens = 'left', x_pred = NULL, miss = NULL,
  tol = 0.0001, M = 10, perc = 0.25, MaxIter = 400, pc = 0.18,
  show_se = TRUE, quiet = FALSE)
}


\arguments{
 \item{cc}{
Vector of censoring indicators of length \code{n}, where \code{n} is the total of observations. For each observation: 0 if non-censored, 1 if censored.
}
\item{x}{
Matrix of covariates of dimension \code{n x l}, where \code{l} is the number of fixed effects including the intercept, if considered (in models which include an intercept \code{x} should contain a column of ones).
}
\item{y}{
Vector of responses of length \code{n}.
}
\item{p}{
Order of the autoregressive process. Must be a positive integer value. For \code{p} equal to \eqn{0} we suggest to use the function \code{\link[SMNCensReg]{CensReg.SMN}} from \code{SMNCensReg} package.
}
\item{cens}{
"left" for left censoring, "right" for right censoring.
}
\item{x_pred}{
Matrix of covariates for responses to be predicted. If x_pred = \code{NULL} no responses are predicted.
}
\item{miss}{
Vector containing the index of missing observations. miss = \code{NULL} indicates that no observations are missing.
}
\item{tol}{
The convergence maximum error permitted.
}
\item{M}{
Size of the Monte Carlo sample generated in each step of the SAEM algorithm. Default=10.
}
\item{perc}{
Percentage of burn-in on the Monte Carlo sample. Default=0.25.
}
\item{MaxIter}{
The maximum number of iterations of the SAEM algorithm. Default=400.
}
\item{pc}{
Percentage of initial iterations of the SAEM algorithm. It is recommended that 50<MaxIter*pc<100. Default=0.18.
}
\item{show_se}{
TRUE or FALSE. Indicates if the standard errors should be estimated. Default=TRUE.
}
\item{quiet}{
TRUE or FALSE. Indicates if printing information should be suppressed. Default=TRUE.
}
}
\details{
The initial values are obtained by ignoring censoring and applying maximum
likelihood estimation with the censored data simply replaced by their censoring limits.
If you want to fit a regression model with autoregressive errors for non-censored data, just set "cc" as a vector of zeros and "cens" as either "right" or "left".
}
\value{
\item{beta}{Estimate of the regression parameters.}
\item{sigma2}{Estimated variance of the white noise process.}
\item{phi}{Estimate of the autoregressive parameters.}
\item{pi1}{Estimate of the first p partial autocorrelations.}
\item{theta}{Vector of parameters estimate (beta, sigma2, phi).}
\item{SE}{Vector of the standard errors of (beta, sigma2, phi).}
\item{loglik}{Log-likelihood value.}
\item{AIC}{Akaike information criterion.}
\item{BIC}{Bayesian information criterion.}
\item{AICcorr}{Corrected Akaike information criterion.}
\item{time}{Processing time.}
\item{pred}{Predicted values (if x_pred is not \code{NULL}).}
\item{criteria}{Attained criteria value.}
\item{yest}{Augmented response variable based on the fitted model.}
\item{yyest}{Final estimative of \code{E(Y\%*\%t(Y))}.}
\item{iter}{Number of iterations until convergence.}
\item{call}{The \code{ARCensReg} call that produced the object.}
\item{tab}{Table of estimates.}
\item{critFin}{Selection criteria.}
\item{cens}{"left" for left censoring, "right" for right censoring.}
\item{nmiss}{Number of missing observations.}
\item{ncens}{Number of censored observations.}
\item{converge}{Logical indicating convergence of the estimation algorithm.}
\item{MaxIter}{The maximum number of iterations used for the SAEM algorithm.}
\item{M}{Size of the Monte Carlo sample generated in each step of the SAEM algorithm.}
\item{pc}{Percentage of burn-in on the Monte Carlo sample.}
\item{timediff}{Time elapsed in processing.}
\item{plot}{A list containing congervence information.}
}

\author{Fernanda L. Schumacher <\email{fernandalschumacher@gmail.com}>,
Victor H. Lachos
<\email{hlachos@ime.unicamp.br}> and Christian E. Galarza <\email{cgalarza88@gmail.com}>

Maintainer: Fernanda L. Schumacher <\email{fernandalschumacher@gmail.com}>
}

\examples{
##simple example (p = l = 1)
#generating a sample
set.seed(23451)
n=50
x=rep(1,n)
dat = rARCens(n=n,beta=2,pit=.5,sig2=.3,x=x,
                cens='left',pcens=.1)

#fitting the model (quick convergence)
fit0 = ARCensReg(dat$data$cc,dat$data$y,x,tol=0.001,
                  pc=.12,M=5,show_se=FALSE)

\dontrun{

##another example (p = l = 2)
#generating a sample
n=100
x=cbind(1,runif(n))
dat = rARCens(n=n,beta=c(2,1),pit=c(.4,-.2),sig2=.5,
                    x=x,cens='left',pcens=.05)

#fitting the model
fit1 = ARCensReg(dat$data$cc,dat$data$y,x,p=2,
                      cens="left",tol=0.0001)
fit1
plot(fit1)

#plotting the augmented variable
par(mfrow=c(1,1))
plot.ts(fit1$yest,lty='dashed',col=4)
lines(dat$data$y)

#simulating missing values
miss = sample(1:n,3)
yMISS = dat$data$y
yMISS[miss] = NA
fit2 = ARCensReg(dat$data$cc,yMISS,x,p=2,miss=miss,
                cens="left",tol=0.0001)
fit2
plot(fit2)

}
}


\references{
Delyon, B., Lavielle, M. & Moulines, E. (1999) Convergence of a stochastic approximation version of the EM algorithm. Journal of the Royal Statistical Society, Series B, 39, 1-38.

Schumacher, F. L., Lachos, V. H. & Dey, D. K. (2016) Censored regression models with autoregressive errors: A likelihood-based perspective. Submitted.
}

\keyword{package}
\keyword{censored regression}
\keyword{autoregressive errors}

\seealso{\code{\link{arima}}, \code{\link[SMNCensReg]{CensReg.SMN}},  \code{\link[ARCensReg]{InfDiag}}}
