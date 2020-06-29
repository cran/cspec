\name{predictiveDFT}
\alias{predictiveDFT}
\title{
Predictive DFT.
}
\description{
Function to calculate the predictive DFT.
}
\usage{
predictiveDFT(x, freq = 2 * pi * (1:length(x))/length(x), taper = FALSE, ar = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
time series vector.
}
  \item{freq}{
frequency vector.
}
  \item{taper}{
logical vector. If TRUE, use the tapered time series to estimate an AR coefficients. Default is FALSE.
}
  \item{ar}{
predetermined AR coefficients to evaluate the predictive DFT. Default is NULL.
}
  \item{\dots}{
Arguments passing ar() function.
}
}
\details{
If ar vector is predetermined, then we evaluate the predictive DFT using the this AR coefficient vector (see eq (2.2) of the reference).
If ar is not specified (default), then we fit the best fitting AR(p) model using AIC and Yule-Walker estimator. 
}
\value{
The complex valued predictive DFT vector.
}
\references{
S. Das, S. Subba Rao, and J. Yang. \emph{Spectral methods for small sample time series: A complete
periodogram approach}. Submitted, 2020.
}
\author{
Junho Yang
}

%\note{ }

\seealso{
\code{\link{complete.pgram}}, \code{\link{tapered.complete.pgram}}
}
\examples{
set.seed(123)
x <- arima.sim(model=list(ar=0.7), n=100)

v1 <- predictiveDFT(x) #default
v2 <- predictiveDFT(x,ar=c(0.7)) #predictive DFT using AR(1) model with coefficient: 0.7.
v3 <- predictiveDFT(x,taper=TRUE) #Using tapered time series to fit the best AR model. 
v4 <- predictiveDFT(x, method="ols") #Using ols method to fit the best AR model.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line