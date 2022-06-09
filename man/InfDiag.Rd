\name{InfDiag}
\alias{InfDiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Influence diagnostic in censored linear regression model with autoregressive errors
}
\description{
It performs influence diagnostic by a local influence approach \insertCite{cook1986assessment}{ARCensReg} with three possible 
perturbation schemes: response perturbation (y), scale matrix perturbation (Sigma), or explanatory 
variable perturbation (x). A benchmark value is calculated that depends on k.
}
\usage{
InfDiag(theta, yest, yyest, x, k = 3, plots = T, 
  indpar = rep(1, length(theta)), perturbation = "y", 
  indcolx = rep(1, ncol(x)))
}
\arguments{
\item{theta}{Vector of estimated parameters.}

\item{yest}{Vector of responses of length \eqn{n} with augmented data. Should be the value yest of 
the \code{\link[ARCensReg]{ARCensReg}} function in the case that at least one observation is censored.}

\item{yyest}{Should be the value yyest of the \code{\link[ARCensReg]{ARCensReg}} function in the case that at least 
one observation is censored. Otherwise, must be \code{y\%*\%t(y)}.}

\item{x}{Matrix of covariates of dimension \eqn{n}x\eqn{l}, where \eqn{l} is the number of fixed effects 
including the intercept, if considered (in models which include an intercept \code{x} should contain a column of ones).}

\item{k}{Constant to be used in the benchmark calculation: \code{M0+k*sd(M0)}.}

\item{plots}{\code{TRUE} or \code{FALSE}. Indicates if a graph should be plotted.}

\item{indpar}{Vector of length equal to the number of parameters, with each element 0 or 1 indicating if the 
respective parameter should be taking into account in the influence calculation.}

\item{perturbation}{Perturbation scheme. Possible values: "y" for response perturbation, "Sigma" for scale matrix 
perturbation or "x" for explanatory variable perturbation.}

\item{indcolx}{If \code{perturbation="x"}, \code{indcolx} must be a vector of length equal to the number of columns of x, 
with each element 0 or 1 indicating if the respective column of x should be perturbed. All columns are perturbed by default.}
}
\details{
The function returns a vector of length \code{n} with the aggregated contribution (\code{M0}) of all eigenvectors 
of the matrix associated with the normal curvature. For details see \insertCite{schumacher2018influence;textual}{ARCensReg}.
}
\value{
\code{M0}
}
\author{
Fernanda L. Schumacher, Katherine L. Valeriano, Victor H. Lachos, Christian E. Galarza,
and Larissa A. Matos
}
\examples{
\donttest{
library(ggplot2)

# Generating the data
set.seed(12341)
x = cbind(1,runif(100))
dat = rARCens(n=100, beta = c(1,-1), pit = c(.4,-.2), sig2=.5, x=x, 
              cens='left', pcens=.05)
# Creating an outlier
dat$data$y[40] = 5
ggplot(dat$data) + geom_line(aes(x=1:100, y=y)) + theme_bw() +
  labs(x="Time")

# Fitting the model
fit = ARCensReg(dat$data$cc, dat$data$lcl, dat$data$ucl, dat$data$y, x, 
                p=2, tol=0.001, show_se=FALSE)

# Influence diagnostic
M0y = InfDiag(theta=fit$res$theta, yest=fit$yest, yyest=fit$yyest,
              x=x, k=3.5, perturbation="y")
M0Sigma = InfDiag(theta=fit$res$theta, yest=fit$yest, yyest=fit$yyest,
                  x=x, k=3.5, perturbation="Sigma")
M0x = InfDiag(theta=fit$res$theta, yest=fit$yest, yyest=fit$yyest,
              x=x, k=3.5, perturbation="x", indcolx=c(0,1))

# Perturbation on a subset of parameters
M0y1 = InfDiag(theta=fit$res$theta, yest=fit$yest, yyest=fit$yyest,
               x=x, k=3.5, perturbation="y", indpar=c(1,1,0,0,0))
M0y2 = InfDiag(theta=fit$res$theta, yest=fit$yest, yyest=fit$yyest,
               x=x, k=3.5, perturbation="y", indpar=c(0,0,1,1,1))
ggplot(data.frame(M0y1,M0y2)) + geom_point(aes(x=M0y1, y=M0y2)) +
  geom_hline(yintercept=mean(M0y2)+3.5*sd(M0y2), linetype="dashed") +
  geom_vline(xintercept=mean(M0y1)+3.5*sd(M0y1), linetype="dashed") +
  theme_bw()}
}
\references{
\insertRef{cook1986assessment}{ARCensReg}

\insertRef{schumacher2018influence}{ARCensReg}

\insertRef{zhu2001local}{ARCensReg}
}
\keyword{package}
\keyword{censored regression}
\keyword{autoregressive errors}
\keyword{influence diagnostic}
\seealso{
\code{\link[ARCensReg]{ARCensReg}}
}
