\name{predict}
\alias{predict.ARpCRM}\alias{predict.ARtpCRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forecast for Autoregressive censored models with Normal and Student-t innovations
}
\description{
Forecast from models fitted by \code{\link[ARCensReg]{ARCensReg}} and \code{\link[ARCensReg]{ARtCensReg}}.
}
\usage{
  \method{predict}{ARpCRM}(object, x_pred, ...)
  
  \method{predict}{ARtpCRM}(object, x_pred, ...)
}
\arguments{
  \item{object}{An object inheriting from class \code{ARpCRM} or \code{ARtpCRM}, representing a fitted AR(p) censored linear model.}
  
  \item{x_pred}{Matrix of covariates for responses to be predicted.} 
  
  \item{...}{Further arguments passed to or from other methods.}
}
\value{
  A time series of predictions.
}
\author{
  Fernanda L. Schumacher, Katherine L. Valeriano, Victor H. Lachos, Christian E. Galarza,
  and Larissa A. Matos
}
\examples{\donttest{
# Generating a sample
set.seed(2839)
n = 210
x = cbind(1, rnorm(n))
dat = rARCens(n=n, beta=c(-1,2), phi=.5, sig2=.3, x=x, cens='left', pcens=.1)

# Fitting the model
data1 = dat$data[1:205,]
fit = ARCensReg(data1$cc, data1$lcl, data1$ucl, data1$y, x[1:205,],
                 M=5, pc=.12, tol=0.001)

# Forecast
y_pred = predict(fit, x[206:n,])
mean((dat$data$y[206:n] - y_pred)^2) # MSPE}
}
\references{
  \insertRef{schumacher2017censored}{ARCensReg}
  
  \insertRef{valeriano2021censored}{ARCensReg}
}
\keyword{package}
\keyword{censored regression}
\keyword{autoregressive errors}
\seealso{
  \code{\link[ARCensReg]{ARCensReg}}, \code{\link[ARCensReg]{ARtCensReg}}
}