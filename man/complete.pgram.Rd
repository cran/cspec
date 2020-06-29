\name{complete.pgram}
\alias{complete.pgram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Complete Periodogram
}
\description{
Function to calculate the complete periodogram.
}
\usage{
complete.pgram(x, freq = 2 * pi * (1:length(x))/length(x), thres=NULL, ...)
}
\arguments{
  \item{x}{
time series vector.
}
  \item{freq}{
frequency vector.
}

  \item{thres}{
(optional) positive threshold value.
}

  \item{\dots}{
Arguments used in the \code{\link{predictiveDFT}} function.
}
}
\details{
The default frequencies are 2*pi*(1:n)/n, where n is a length of a time series.
}
\value{
Real part of the complete periodogram vector. 
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
\code{\link{predictiveDFT}}
}
\examples{
set.seed(123)
x <- arima.sim(model=list(ar=0.7), n=100)
v <- complete.pgram(x)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line