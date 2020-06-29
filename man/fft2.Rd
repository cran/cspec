\name{fft2}
\alias{fft2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A new DFT function
}
\description{
Function to calculate the DFT of time series vector as in the reference.
}
\usage{
fft2(x, freq = 2 * pi * (1:length(x))/length(x))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
time series vector.
}
  \item{freq}{
frequency vector.
}
}
\details{
sqrt(n) standardized. The default frequencies are 2*pi*(1:n)/n, where n is a length of a time series.
}
\value{
DFT vector.
}
\references{
S. Das, S. Subba Rao, and J. Yang. \emph{Spectral methods for small sample time series: A complete
periodogram approach}. Submitted, 2020.
}
\author{
Junho Yang
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{taperDFT}}
}
\examples{
set.seed(123)
x <- arima.sim(model=list(ar=0.7), n=100)
v <- fft2(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line