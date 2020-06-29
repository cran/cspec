\name{completeDFT}
\alias{completeDFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Complete DFT
}
\description{
Function to calculate the complete DFT.
}
\usage{
completeDFT(x, freq = 2 * pi * (1:length(x))/length(x), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
time series vector.
}
  \item{freq}{
frequency vector.
}
  \item{\dots}{
Arguments used in the \code{\link{predictiveDFT}} function.
}
}
\details{
The default frequencies are 2*pi*(1:n)/n, where n is a length of a time series.
}
\value{
Complex valued complete DFT vector. 
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


\seealso{
\code{\link{predictiveDFT}}
}
\examples{
set.seed(123)
x <- arima.sim(model=list(ar=0.7), n=100)
v <- completeDFT(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line