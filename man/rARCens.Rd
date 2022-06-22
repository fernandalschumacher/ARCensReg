\name{rARCens}
\alias{rARCens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating censored autoregressive data
}
\description{
It simulates a censored response variable with autoregressive errors of order \eqn{p} following normal or Student-t innovations, with an established censoring rate.
}
\usage{
rARCens(n, beta, phi, sig2 = 1, x = rep(1, n), cens = "left", 
  pcens = 0.1, innov = "norm", nu = NULL)
}
\arguments{
\item{n}{Length of the desired time serie.}

\item{beta}{Vector of theoretical regression parameters of length \eqn{l}.}

\item{phi}{Vector of theoretical autoregressive coefficients of length \eqn{p}.}

\item{sig2}{Theoretical variance of the error.}

\item{x}{Matrix of covariates of dimension \eqn{n}x\eqn{l} (in models that include 
an intercept \code{x} should contain a column of ones).}

\item{cens}{\code{'left'} for left censoring, \code{'right'} for right censoring.}

\item{pcens}{Desired censoring rate.}

\item{innov}{Distribution of the innovation variable. The values are \code{'norm'} and \code{'t'} for normal and Student-t distribution, respectively.}

\item{nu}{Degrees of freedom for Student-t innovations.}
}
\value{
\item{data}{Generated response (y), censoring indicator (cc), and lower (lcl) and upper (ucl) bounds of the interval, which contains the true value of the censored observation.}
\item{param}{Theoretical parameters (beta, sig2, phi).}
}
\note{
For data generation with Student-t innovations, the first \eqn{p} observations are not censored.
}
\author{
Fernanda L. Schumacher, Katherine L. Valeriano, Victor H. Lachos, Christian E. Galarza,
and Larissa A. Matos
}
\examples{
library(ggplot2)

## Example 1: Generating a sample with normal innovations
set.seed(1234)
dat = rARCens(n=100, beta=c(1,-1), phi=c(.48,-.2), sig2=.5,
              x=cbind(1,runif(100)), cens='left', pcens=.10)

# Plotting the time serie
ggplot(data.frame(dat$data$y), aes(x=1:100, y=dat$data$y)) + geom_line() + 
  geom_line(aes(x=1:100, y=dat$data$ucl), color="red", linetype="twodash") + 
  labs(x="Time", y=bquote(y["obs"])) + theme_bw()

table(dat$data$cc)

dat$param
#[1]  1.00 -1.00  0.50  0.48 -0.20

## Example 2: Generating a sample with Student-t innovations
set.seed(8278)
dat1 = rARCens(n=100, beta=c(1,-1), phi=c(.48,-.2), sig2=.5,
               x=cbind(1,rnorm(100)), cens='right', pcens=.10, 
               innov='t', nu=3)

# Plotting the time serie
ggplot(data.frame(dat1$data$y), aes(x=1:100, y=dat1$data$y)) + geom_line() + 
  geom_line(aes(x=1:100, y=dat1$data$lcl), color="red", linetype="twodash") + 
  labs(x="Time", y=bquote(y["obs"])) + theme_bw()
  
dat1$param
#[1]  1.00 -1.00  0.50  0.48 -0.20  3.00
}
\keyword{package}
\keyword{censored regression}
\keyword{autoregressive errors}
\seealso{
\code{\link[ARCensReg]{ARCensReg}}, \code{\link[ARCensReg]{ARtCensReg}}
}
