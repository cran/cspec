\name{tapered.complete.pgram}
\alias{tapered.complete.pgram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tapered Complete Periodogram
}
\description{
Function to calculate the tapered complete periodogram.
}
\usage{
tapered.complete.pgram(x, freq=2*pi*(1:length(x))/length(x), taperx=NULL, thres=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
time series vector.
}
  \item{freq}{
frequency vector.
}
  \item{taperx}{
predetermined tapered DFT. If NULL, we use default settings of the taperedDFT function to calculate the tapered DFT.
}

  \item{thres}{
positive threshold value.
}

  \item{\dots}{
Arguments used in the \code{\link{predictiveDFT}} function.
}
}
\details{
The default frequencies are 2*pi*(1:n)/n, where n is a length of a time series.
}
\value{
Real part of the tapered complete periodogram vector.
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
\code{\link{predictiveDFT}}, \code{\link{complete.pgram}}
}
\examples{
set.seed(123)
x <- arima.sim(model=list(ar=0.7), n=100)
v <- tapered.complete.pgram(x)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line