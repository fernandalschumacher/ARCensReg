\docType{data}
\name{phosphorus}
\alias{phosphorus}
\title{
Phosphorus concentration data
}
\description{
The phosphorus concentration (P) data of West Fork Cedar River at Finchford, 
Iowa, USA, collected under the ambient water quality program conducted by the
Iowa Department of Natural Resources (Iowa DNR), were observed monthly from 
10/1998 to 10/2013 (n=181). The phosphorus concentration measurement was subject
to a detection limit (lcl); thereby, the P data are left-censored. The dataset
was first available in the R package carx. 

The water discharge dataset was obtained from the website of the U.S. Geological Survey (site number 05458900), and it is measured in cubic feet per second.
}
\usage{
data(phosphorus)
}
\format{
This data frame contains the following columns:
\describe{
  \item{\code{lP}}{Logarithm of the phosphorus concentration.}
  \item{\code{cc}}{Left censoring indicator (1 if the observation is left-censored and 0 otherwise).}
  \item{\code{lQ}}{Logarithm of the water discharge.}
  \item{\code{lcl}}{Censoring limit.}
  \item{\code{time}}{Year-Month.}
}
}
\source{
\url{https://waterdata.usgs.gov/ia/nwis/monthly/}

\url{https://CRAN.R-project.org/package=carx}
}
\examples{\donttest{
library(ggplot2)

data(phosphorus)
n = nrow(phosphorus)

ggplot(phosphorus) + geom_line(aes(x=1:n, y=lP)) + 
  geom_line(aes(x=1:n, y=lcl), color="red", linetype="dashed") + 
  labs(x="Time") + theme_bw()

# Proportion of censoring
prop.table(table(phosphorus$cc))

# A censored regression model
x   = cbind(1, phosphorus$lQ)
cc  = phosphorus$cc
lcl = rep(-Inf, n)
ucl = phosphorus$lcl
miss =  which(is.na(phosphorus$lP))
cc[miss]  = 1 
ucl[miss] = Inf

# Fitting a model with normal innovations
set.seed(8765)
mod1 = ARCensReg(cc, lcl, ucl, phosphorus$lP, x, p=1, tol=.001)

# Fitting a model with Student-t innovations
set.seed(287399)
mod2 = ARtCensReg(cc, lcl, ucl, phosphorus$lP, x, p=1, tol=.001)

# Plotting observed and imputed values
data.plot = data.frame(y=phosphorus$lP, ynorm=mod1$yest, yt=mod2$yest)
#
ggplot(data.plot) + geom_line(aes(x=1:n, y=ynorm), color=4) +
  geom_line(aes(x=1:n, y=yt), color="deeppink", linetype="dashed") + 
  geom_line(aes(x=1:n, y=y)) + labs(x="Time", y="lP") + theme_bw()

# Imputed values
data.plot[cc==1,]
}}
\seealso{
\code{\link[ARCensReg]{ARCensReg}}, \code{\link[ARCensReg]{ARtCensReg}}
}
\keyword{datasets}
