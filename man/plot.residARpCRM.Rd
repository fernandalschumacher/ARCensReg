\name{plot.residARpCRM}
\alias{plot.residARpCRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Show diagnostic residual plots
}
\description{
  It returns four plots for the quantile residuals: the time series plot of the residuals, the quantile-quantile plot, the histogram, and the ACF plot of the residuals.
}

\usage{
  \method{plot}{residARpCRM}(x, \dots)
}

\arguments{
  \item{x}{An object inheriting from class \code{residARpCRM} obtained as an output of function \code{\link{residuals}}.}
  \item{\dots}{Additional arguments.}
}

\value{
  A ggplot object.
}

\author{
Fernanda L. Schumacher, Katherine L. Valeriano, Victor H. Lachos, Christian E. Galarza, and Larissa A. Matos
}


\seealso{
  \code{\link[ggplot2]{ggplot}}, \code{\link{ARCensReg}}, \code{\link{ARtCensReg}}, \code{\link{residuals.ARpCRM}}, \code{\link{residuals.ARtpCRM}}
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
r1 = residuals(mod1)
class(r1)
plot(r1)

# Fitting the model with Student-t innovations
mod2 = ARtCensReg(dat1$data$cc, dat1$data$lcl, dat1$data$ucl, dat1$data$y, 
                  x, p=2, tol=0.001)
r2 = residuals(mod2)
plot(r2)


## Example 2: Generating heavy-tailed data
set.seed(12341)
x = cbind(1, runif(300))
dat2 = rARCens(n=300, beta=c(1,-1), phi=c(.48,-.2), sig2=.5, x=x, 
              cens='left', pcens=.05, innov="t", nu=3)

# Fitting the model with normal innovations
mod3 = ARCensReg(dat2$data$cc, dat2$data$lcl, dat2$data$ucl, dat2$data$y,
                 x, p=2, tol=0.001)
r3 = residuals(mod3)
plot(r3)

# Fitting the model with Student-t innovations
mod4 = ARtCensReg(dat2$data$cc, dat2$data$lcl, dat2$data$ucl, dat2$data$y,
                  x, p=2, tol=0.001)
r4 = residuals(mod4)
plot(r4)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
