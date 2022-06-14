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
InfDiag(object, k = 3, indpar = rep(1, length(object$theta)), 
  indcolx = rep(1, ncol(object$x)), perturbation = "y")
}
\arguments{
\item{object}{Object of class \code{'ARpCRM'} given as an output of function \code{\link[ARCensReg]{ARCensReg}}.}

\item{k}{Constant to be used in the benchmark calculation: \code{M0+k*sd(M0)}.}

\item{indpar}{Vector of length equal to the number of parameters, with each element 0 or 1 indicating if the 
respective parameter should be considered in the influence calculation.}

\item{indcolx}{If \code{perturbation="x"}, \code{indcolx} must be a vector of length equal to the number of columns of x, 
with each element 0 or 1 indicating if the respective column of x should be perturbed. All columns are perturbed by default.}

\item{perturbation}{Perturbation scheme. Possible values: "y" for response perturbation, "Sigma" for scale matrix 
perturbation, or "x" for explanatory variable perturbation.}
}
\details{
The function returns a vector of length \code{n} with the aggregated contribution (\code{M0}) of all eigenvectors 
of the matrix associated with the normal curvature. For details see \insertCite{schumacher2018influence;textual}{ARCensReg}.
}
\value{
An object of class "DiagARpCRM" with the following components is returned:
\item{M0}{Vector of length \eqn{n} with the aggregated contribution of all eigenvectors of the matrix associated with the normal curvature.}
\item{perturbation}{Perturbation scheme.}
\item{benchmark}{\code{M0 + k*sd(M0)}.}
}
\author{
Fernanda L. Schumacher, Katherine L. Valeriano, Victor H. Lachos, Christian E. Galarza,
and Larissa A. Matos
}
\examples{\donttest{
library(ggplot2)

# Generating the data
set.seed(12341)
x = cbind(1,runif(100))
dat = rARCens(n=100, beta=c(1,-1), phi=c(.48,-.2), sig2=.5, x=x, 
              cens='left', pcens=.05)
              
# Creating an outlier
dat$data$y[40] = 5
ggplot(dat$data) + geom_line(aes(x=1:100, y=y)) + theme_bw() +
  labs(x="Time")

# Fitting the model
fit = ARCensReg(dat$data$cc, dat$data$lcl, dat$data$ucl, dat$data$y, x, 
                p=2, tol=0.001, show_se=FALSE)

# Influence diagnostic
M0y = InfDiag(fit, k=3.5, perturbation="y")
plot(M0y)
M0Sigma = InfDiag(fit, k=3.5, perturbation="Sigma")
plot(M0Sigma)
M0x = InfDiag(fit, k=3.5, indcolx=c(0,1), perturbation="x")
plot(M0x)

# Perturbation on a subset of parameters
M0y1 = InfDiag(fit, k=3.5, indpar=c(1,1,0,0,0), perturbation="y")$M0
M0y2 = InfDiag(fit, k=3.5, indpar=c(0,0,1,1,1), perturbation="y")$M0
#
ggplot(data.frame(M0y1,M0y2)) + geom_point(aes(x=M0y1, y=M0y2)) +
  geom_hline(yintercept=mean(M0y2)+3.5*sd(M0y2), linetype="dashed") +
  geom_vline(xintercept=mean(M0y1)+3.5*sd(M0y1), linetype="dashed") +
  theme_bw()
}}
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
