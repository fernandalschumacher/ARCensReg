\name{residuals}
\alias{residuals.ARpCRM}\alias{residuals.ARtpCRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract model residuals from ARpCRM and ARtpCRM objects
}
\description{
The conditional residuals are obtained by subtracting the fitted values from the response vector, while the quantile residuals are obtained by inverting the estimated distribution function for each observation to obtain approximately normally distributed residuals. See, for instance, \insertCite{dunn1996randomized;textual}{ARCensReg} and \insertCite{kalliovirta2012misspecification;textual}{ARCensReg}.
}
\usage{
\method{residuals}{ARpCRM}(object, ...)

\method{residuals}{ARtpCRM}(object, ...)
}
\arguments{
\item{object}{An object inheriting from class \code{ARpCRM} or \code{ARtpCRM}, representing a fitted AR(p) censored linear model.}
\item{...}{Further arguments passed to or from other methods.}
}
\value{
An object of class "residARpCRM", with the following components:
\item{residuals}{Vector with the conditional residuals of length \eqn{n}.}
\item{quantile.resid}{Vector with the quantile residuals of length \eqn{n}.}
Generic function plot has methods to show a dispersion plot (residual vs. time), an autocorrelation plot, a histogram, and Quantile-Quantile (Q-Q) plot for the quantile residuals.
}
\author{
Fernanda L. Schumacher, Katherine L. Valeriano, Victor H. Lachos, Christian E. Galarza,
and Larissa A. Matos
}
\examples{\donttest{
## Example 1: Generating data with normal innovations
set.seed(93899)
x = cbind(1, runif(300))
dat1 = rARCens(n=300, beta=c(1,-1), phi=c(.48,-.2), sig2=.5, x=x, 
              cens='left', pcens=.05, innov="norm")

# Fitting the model with normal innovations
mod1 = ARCensReg(dat1$data$cc, dat1$data$lcl, dat1$data$ucl, dat1$data$y, 
                 x, p=2, tol=0.001)
mod1$tab
plot(residuals(mod1))

# Fitting the model with Student-t innovations
mod2 = ARtCensReg(dat1$data$cc, dat1$data$lcl, dat1$data$ucl, dat1$data$y, 
                  x, p=2, tol=0.001)
mod2$tab
plot(residuals(mod2))


## Example 2: Generating heavy-tailed data
set.seed(12341)
x = cbind(1, runif(300))
dat2 = rARCens(n=300, beta=c(1,-1), phi=c(.48,-.2), sig2=.5, x=x, 
              cens='left', pcens=.05, innov="t", nu=3)

# Fitting the model with normal innovations
mod3 = ARCensReg(dat2$data$cc, dat2$data$lcl, dat2$data$ucl, dat2$data$y,
                 x, p=2, tol=0.001)
mod3$tab
plot(residuals(mod3))

# Fitting the model with Student-t innovations
mod4 = ARtCensReg(dat2$data$cc, dat2$data$lcl, dat2$data$ucl, dat2$data$y,
                  x, p=2, tol=0.001)
mod4$tab
plot(residuals(mod4))
}}
\references{
\insertRef{dunn1996randomized}{ARCensReg}

\insertRef{kalliovirta2012misspecification}{ARCensReg}
}
\keyword{package}
\keyword{censored regression}
\keyword{autoregressive errors}
\seealso{
\code{\link[ARCensReg]{ARCensReg}}, \code{\link[ARCensReg]{ARtCensReg}}
}