\docType{data}
\name{CloudCeiling}
\alias{CloudCeiling}
\title{
Cloud ceiling height
}
\description{
The cloud ceiling heights, collected by the National Center for Atmospheric Research (NCAR), 
were observed hourly in San Francisco during March 1989, consisting of n=716 observations 
\insertCite{park2007censored}{ARCensReg}.
}
\usage{
data(CloudCeiling)
}
\format{
This data frame contains the following columns:
\describe{
    \item{\code{y}}{Logarithm of the cloud ceiling heights.}
    \item{\code{cc}}{Right censoring indicator (1 if the observation is right-censored and 0 otherwise).}
}
}
\source{
\insertRef{park2007censored}{ARCensReg}
}
\examples{
library(ggplot2)

data(CloudCeiling)
ggplot(CloudCeiling) + geom_line(aes(x=1:length(y), y=y)) + 
  labs(x="Time") + theme_bw()

# Proportion of censoring
prop.table(table(CloudCeiling$cc))

\dontrun{
# A censored regression model
## This may take a long time due to the number of censored observations.
## For other examples see help(ARCensReg).

x   = as.matrix(rep(1, length(CloudCeiling$y)))
cc  = CloudCeiling$cc
lcl = CloudCeiling$y
ucl = rep(Inf, length(CloudCeiling$y))
miss =  which(is.na(CloudCeiling$y))
cc[miss]  = 1
lcl[miss] = -Inf
AR_reg = ARCensReg(cc, lcl, ucl, CloudCeiling$y, x, p=1, tol=.001)}
}
\seealso{
\code{\link{ARCensReg}}, \code{\link{ARtCensReg}}
}
\keyword{datasets}
