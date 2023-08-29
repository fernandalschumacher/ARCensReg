\name{ARtCensReg}
\alias{ARtCensReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Censored autoregressive regression model with Student-t innovations
}
\description{
It fits a univariate left, right, or interval censored linear regression model 
with autoregressive errors considering Student-t innovations, through the SAEM algorithm. 
It provides estimates and standard errors of the parameters, supporting missing 
values on the dependent variable. 
}
\usage{
ARtCensReg(cc, lcl = NULL, ucl = NULL, y, x, p = 1, M = 10,
  perc = 0.25, MaxIter = 400, pc = 0.18, nufix = NULL, tol = 1e-04,
  show_se = TRUE, quiet = FALSE)
}
\arguments{
\item{cc}{Vector of censoring indicators of length \eqn{n}, where \eqn{n} is the 
total observations. For each observation: 0 if non-censored, 1 if censored/missing.}


\item{lcl, ucl}{Vectors of length \eqn{n} that represent the lower and upper bounds 
of the interval, which contains the observade value of the censored observation. 
Default=\code{NULL}, indicating no-censored data. See details for more information.}

\item{y}{Vector of responses of length \eqn{n}.}
  
\item{x}{Matrix of covariates of dimension \eqn{n \times l}, where \eqn{l} is the 
number of fixed effects including the intercept, if considered (in models which 
include an intercept, \code{x} should contain a column of ones).}
  
\item{p}{Order of the autoregressive process. It must be a positive integer value.}
  
\item{M}{Size of the Monte Carlo sample generated in each step of the SAEM algorithm.
Default=10.}
  
\item{perc}{Percentage of burn-in on the Monte Carlo sample. Default=0.25.}
  
\item{MaxIter}{The maximum number of iterations of the SAEM algorithm. Default=400.}
  
\item{pc}{Percentage of initial iterations of the SAEM algorithm with no memory. It is recommended that 
\code{50<MaxIter*pc<100}. Default=0.18.}

\item{nufix}{If the degrees of freedom (\eqn{\nu}) are unknown, \code{nufix} 
should be equal to \code{NULL}; otherwise, it must be a number greater than 2.}

\item{tol}{The convergence maximum error permitted.}
  
\item{show_se}{\code{TRUE} or \code{FALSE}. Indicates if the standard errors should be 
estimated. Default=\code{TRUE}.}
  
\item{quiet}{\code{TRUE} or \code{FALSE}. Indicates if printing information should be 
suppressed. Default=\code{FALSE}.}
}
\details{
The linear regression model with autocorrelated errors, defined as a discrete-time autoregressive (AR) process of order \eqn{p}, at time \eqn{t} is given by

\eqn{Y_t = x_t^T \beta + \xi_t,}

\eqn{\xi_t = \phi_1 \xi_{t-1} + ... + \phi_p \xi_{t-p} + \eta_t, t=1,..., n,}

where \eqn{Y_t} is the response variable, \eqn{\beta = (\beta_1,..., \beta_l)^T} is 
a vector of regression parameters of dimension \eqn{l}, \eqn{x_t = (x_{t1},..., x_{tl})^T} 
is a vector of non-stochastic regressor variables values, and \eqn{\xi_t} is the AR 
error with \eqn{\eta_t} being a shock of disturbance following the Student-t distribution 
with \eqn{\nu} degrees of freedom, \eqn{\phi = (\phi_1,..., \phi_p)^T} being the vector of AR 
coefficients, and \eqn{n} denoting the sample size. 

It is assumed that \eqn{Y_t} is not fully observed for all \eqn{t}. 
For left censored observations, we have \code{lcl=-Inf} and \code{ucl=}\eqn{V_t}, 
such that the true value \eqn{Y_t \leq V_t}. For right censoring, \code{lcl=}\eqn{V_t} 
and \code{ucl=Inf}, such that \eqn{Y_t \geq V_t}. For interval censoring, \code{lcl} and \code{ucl} must be finite values, such that \eqn{V_{1t} \leq Y_t \leq V_{2t}}. Missing data can be defined by setting \code{lcl=-Inf} and \code{ucl=Inf}.

The initial values are obtained by ignoring censoring and applying maximum likelihood 
estimation with the censored data replaced by their censoring limits. Moreover, 
just set \code{cc} as a vector of zeros to fit a regression model with autoregressive 
errors for non-censored data.
}
\section{Warning }{This algorithm assumes that the first \eqn{p} values in the response vector are completely observed.
}
\value{
An object of class "ARtpCRM" representing the AR(p) censored regression Student-t fit. Generic functions such as print and summary have methods to show the results of the fit. The function plot provides convergence graphics for the parameter estimates.

Specifically, the following components are returned:
\item{beta}{Estimate of the regression parameters.}
\item{sigma2}{Estimated scale parameter of the innovation.}
\item{phi}{Estimate of the autoregressive parameters.}
\item{nu}{Estimated degrees of freedom.}
\item{theta}{Vector of parameters estimate (\eqn{\beta, \sigma^2, \phi, \nu}).}
\item{SE}{Vector of the standard errors of (\eqn{\beta, \sigma^2, \phi, \nu}).}
\item{yest}{Augmented response variable based on the fitted model.}
\item{uest}{Final estimated weight variables.}
\item{x}{Matrix of covariates of dimension \eqn{n \times l}.}
\item{iter}{Number of iterations until convergence.}
\item{criteria}{Attained criteria value.}
\item{call}{The \code{ARtCensReg} call that produced the object.}
\item{tab}{Table of estimates.}
\item{cens}{"left", "right", or "interval" for left, right, or interval censoring, respectively.}
\item{nmiss}{Number of missing observations.}
\item{ncens}{Number of censored observations.}
\item{converge}{Logical indicating convergence of the estimation algorithm.}
\item{MaxIter}{The maximum number of iterations used for the SAEM algorithm.}
\item{M}{Size of the Monte Carlo sample generated in each step of the SAEM algorithm.}
\item{pc}{Percentage of initial iterations of the SAEM algorithm with no memory.}
\item{time}{Time elapsed in processing.}
\item{plot}{A list containing convergence information.}
}
\author{
Katherine L. Valeriano, Fernanda L. Schumacher, and Larissa A. Matos
}
\examples{
## Example 1: (p = l = 1)
# Generating a sample
set.seed(1234)
n = 80
x = rep(1, n)
dat = rARCens(n=n, beta=2, phi=.6, sig2=.3, x=x, cens='right', pcens=.05, 
              innov='t', nu=4)

# Fitting the model (quick convergence)
fit0 = ARtCensReg(dat$data$cc, dat$data$lcl, dat$data$ucl, dat$data$y, x,
                 M=5, pc=.12, tol=0.001)
fit0
\donttest{
## Example 2: (p = l = 2)
# Generating a sample
set.seed(783796)
n = 200
x = cbind(1, runif(n))
dat = rARCens(n=n, beta=c(2,1), phi=c(.48,-.2), sig2=.5, x=x, cens='left',
              pcens=.05, innov='t', nu=5)
  
# Fitting the model with nu known
fit1 = ARtCensReg(dat$data$cc, dat$data$lcl, dat$data$ucl, dat$data$y, x,
                  p=2, M=15, pc=.20, nufix=5)
summary(fit1)
plot(fit1)

# Fitting the model with nu unknown
fit2 = ARtCensReg(dat$data$cc, dat$data$lcl, dat$data$ucl, dat$data$y, x,
                 p=2, M=15, pc=.20)
summary(fit2)
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
\code{\link{arima}}, \code{\link[ARCensReg]{ARCensReg}}
}