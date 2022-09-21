\name{plot.DiagARpCRM}
\alias{plot.DiagARpCRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot influence diagnostic measures
}
\description{
  Plot method for objects of class "DiagARpCRM".
}

\usage{
  \method{plot}{DiagARpCRM}(x, \dots)
}

\arguments{
  \item{x}{An object inheriting from class \code{DiagARpCRM}. The influence diagnostic measures are calculated by function \code{\link{InfDiag}}, with three possible perturbation schemes: response perturbation (y), scale matrix perturbation (Sigma), or explanatory variable perturbation (x).}
  \item{\dots}{Additional arguments.}
}

\value{
  A ggplot object, plotting the index versus the influence diagnostic measure.
}

\author{
Fernanda L. Schumacher, Katherine L. Valeriano, Victor H. Lachos, Christian E. Galarza, and Larissa A. Matos
}


\seealso{
  \code{\link{ggplot}}, \code{\link{InfDiag}}, \code{\link{ARCensReg}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")