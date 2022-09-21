\name{plot}
\alias{plot.ARpCRM} \alias{plot.ARtpCRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot an ARpCRM or ARtpCRM object
  %%  ~~function to do ... ~~
}
\description{
It displays convergence graphs for the parameters estimates (for the case with at 
least one censored observation). The dashed line indicates the iteration of the 
SAEM algorithm that simulations start being smoothed. 
}

\usage{
  \method{plot}{ARpCRM}(x, \dots)
  
  \method{plot}{ARtpCRM}(x, \dots)
}

\arguments{
  \item{x}{An object inheriting from class \code{ARpCRM} or \code{ARtpCRM}, 
  representing a fitted censored autoregressive model of order \eqn{p}, with 
  normal and Student-t innovations, respectively.}
  \item{\dots}{Additional arguments.}
}

\value{
  A ggplot object.
}

\author{
  Fernanda L. Schumacher, Katherine L. Valeriano, Victor H. Lachos, Christian E. Galarza, and Larissa A. Matos
}


\seealso{
\code{\link{ggplot}}, \code{\link{ARCensReg}}, \code{\link{ARtCensReg}}
}
\examples{
n = 50; x = rep(1, n)
dat = rARCens(n=n, beta=2, phi=.5, sig2=.3, x=x, cens='left', pcens=.1)

fit = ARCensReg(dat$data$cc, dat$data$lcl, dat$data$ucl, dat$data$y, x,
                M=5, pc=.12, tol=0.001, show_se=FALSE)
plot(fit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")