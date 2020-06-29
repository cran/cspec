\name{taperDFT}
\alias{taperDFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tapered DFT
}
\description{
Function to calculate the tapered DFT.
}
\usage{
taperDFT(x, freq = 2 * pi * (1:length(x))/length(x), regularization.type = "1", p = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
time series vector.
}
  \item{freq}{
frequency vector.
}
  \item{regularization.type}{
character. if "1": sum of the taper equals to n, if "2": square sum of the taper equals to n. Default is "1".
}
  \item{p}{
proprotion of taper.
}
}
\details{
Tapered DFT. The default taper is a Tukey's (or Cosine-bell) taper. 
}
\value{
Complex valued tapered DFT vector.
}

%\references{
%}

\author{
Junho Yang.
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fft2}}
}
\examples{
set.seed(123)
x <- arima.sim(model=list(ar=0.7), n=100)
v <- taperDFT(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
