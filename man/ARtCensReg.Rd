\name{ARtCensReg}
\alias{ARtCensReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Censored autoregressive regression model with Student-t innovations
}
\description{
It fits a univariate left, right, or interval censored linear regression model 
with autoregressive errors considering Student-t innovations through the SAEM algorithm. 
It provides estimates and standard errors of the parameters, predicts future 
observations, and supports missing values on the dependent variable. 
}
\usage{
ARtCensReg(cc, lcl = NULL, ucl = NULL, y, x, p = 1, x_pred = NULL, 
  tol = 1e-04, M = 10, perc = 0.25, MaxIter = 400, pc = 0.18, 
  nufix = NULL, show_se = TRUE, quiet = FALSE)
}
\arguments{
\item{cc}{Vector of censoring indicators of length \eqn{n}, where \eqn{n} is the 
total of observations. For each observation: 0 if non-censored, 1 if censored/missing.}
  
\item{lcl, ucl}{Vectors of length \eqn{n}, representing the lower and upper bounds of the interval which contains the true value of the censored observation. Default =\code{NULL}, indicating no-censored data.}
  
\item{y}{Vector of responses of length \eqn{n}.}
  
\item{x}{Matrix of covariates of dimension \eqn{n}x\eqn{l}, where \eqn{l} is the number 
of fixed effects including the intercept, if considered (in models which include an 
intercept \code{x} should contain a column of ones).}
  
\item{p}{Order of the autoregressive process. Must be a positive integer value.}
  
\item{x_pred}{Matrix of covariates for responses to be predicted. If x_pred = \code{NULL} 
no responses are predicted.}
  
\item{tol}{The convergence maximum error permitted.}
  
\item{M}{Size of the Monte Carlo sample generated in each step of the SAEM algorithm.
Default=10.}
  
\item{perc}{Percentage of burn-in on the Monte Carlo sample. Default=0.25.}
  
\item{MaxIter}{The maximum number of iterations of the SAEM algorithm. Default=400.}
  
\item{pc}{Percentage of initial iterations of the SAEM algorithm with no-memory. It is recommended that 
\code{50<MaxIter*pc<100}. Default=0.18.}

\item{nufix}{If the degrees of freedom (\eqn{\nu}) are unknown \code{nufix} should be equal \code{NULL}, otherwise 
it must be a number greater than 2.}
  
\item{show_se}{\code{TRUE} or \code{FALSE}. Indicates if the standard errors should be 
estimated. Default=\code{TRUE}.}
  
\item{quiet}{\code{TRUE} or \code{FALSE}. Indicates if printing information should be 
suppressed. Default=\code{FALSE}.}
}
\details{
The linear regression model with autocorrelated errors, defined as a discrete-time autoregressive (AR) process of order \eqn{p}, at time \eqn{t} is given by

\eqn{Y_t = x_t^T \beta + \xi_t,}

\eqn{\xi_t = \phi_1 \xi_{t-1} + ... + \phi_p \xi_{t-p} + \eta_t, t=1,..., n,}

where \eqn{Y_t} is the response variable (we assume that \eqn{Y_t} is not fully observed for all \eqn{t}), \eqn{\beta = (\beta_1,..., \beta_l)^T} is a vector of regression parameters of dimension \eqn{l}, \eqn{x_t = (x_{t1},..., x_{tl})^T} is a vector of non-stochastic regressor variables values, and \eqn{\xi_t} is the AR error with \eqn{\eta_t} being a shock of disturbance following the Student-t distribution, \eqn{\phi = (\phi_1,..., \phi_p)^T} being the vector of AR coefficients, and \eqn{n} denoting the sample size. 

The initial values are obtained by ignoring censoring and applying maximum likelihood 
estimation with the censored data replaced by their censoring limits. To fit a 
regression model with autoregressive errors for non-censored data, just set
\code{cc} as a vector of zeros.
}
\section{Warning }{This algorithm assumes that the first \eqn{p} values in the response vector are completely observed.
}
\value{
\item{beta}{Estimate of the regression parameters.}
\item{sigma2}{Estimated scale parameter of the innovation.}
\item{phi}{Estimate of the autoregressive parameters.}
\item{nu}{Estimated degrees of freedom.}
\item{theta}{Vector of parameters estimate (\eqn{\beta, \sigma^2, \phi, \nu}).}
\item{SE}{Vector of the standard errors of (\eqn{\beta, \sigma^2, \phi, \nu}).}
\item{pred}{Predicted values (if x_pred is not \code{NULL}).}
\item{criteria}{Attained criteria value.}
\item{yest}{Augmented response variable based on the fitted model.}
\item{uest}{Final estimated weight variables.}
\item{x}{Matrix of covariates of dimension \eqn{n}x\eqn{l}.}
\item{iter}{Number of iterations until convergence.}
\item{call}{The \code{ARtCensReg} call that produced the object.}
\item{tab}{Table of estimates.}
\item{cens}{"left", "right", or "interval" for left, right, or interval censoring, respectively.}
\item{nmiss}{Number of missing observations.}
\item{ncens}{Number of censored observations.}
\item{converge}{Logical indicating convergence of the estimation algorithm.}
\item{MaxIter}{The maximum number of iterations used for the SAEM algorithm.}
\item{M}{Size of the Monte Carlo sample generated in each step of the SAEM algorithm.}
\item{pc}{Percentage of initial iterations of the SAEM algorithm with no-memory.}
\item{time}{Time elapsed in processing.}
\item{plot}{A list containing congervence information.}
}
\author{
Fernanda L. Schumacher, Katherine L. Valeriano, Victor H. Lachos, Christian E. Galarza, and Larissa A. Matos
}
\examples{
## Simple example (p = l = 1)
# Generating a sample
set.seed(1234)
n = 80
x = rep(1, n)
dat = rARCens(n=n, beta=2, pit=.6, sig2=.3, x=x, cens='right', pcens=.05, 
              innov='t', nu=4)

# Fitting the model (quick convergence)
fit0 = ARtCensReg(dat$data$cc, dat$data$lcl, dat$data$ucl, dat$data$y, x,
                 tol=0.001, M=5, pc=.12)
\donttest{
## Another example (p = l = 2)
# Generating a sample
set.seed(783796)
n = 200
x = cbind(1, runif(n))
dat = rARCens(n=n, beta=c(2,1), pit=c(.4,-.2), sig2=.5, x=x, cens='left',
              pcens=.05, innov='t', nu=5)
  
# Fitting the model with nu known
fit1 = ARtCensReg(dat$data$cc, dat$data$lcl, dat$data$ucl, dat$data$y, x,
                  p=2, M=15, pc=.20, nufix=5)
fit1$tab
plot(fit1)

# Fitting the model with nu unknown
fit2 = ARtCensReg(dat$data$cc, dat$data$lcl, dat$data$ucl, dat$data$y, x,
                 p=2, M=15, pc=.20)
fit2$tab
plot(fit2)}
}
\references{
  \insertRef{delyon1999convergence}{ARCensReg}
  
  \insertRef{valeriano2021censored}{ARCensReg}
}
\keyword{package}
\keyword{censored regression}
\keyword{autoregressive errors}
\seealso{
\code{\link{arima}}, \code{\link[ARCensReg]{ARCensReg}}, 
\code{\link[SMNCensReg]{CensReg.SMN}}
}