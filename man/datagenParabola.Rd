\name{datagenParabola}
\alias{datagenParabola}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Parabola generating function
}
\description{
Newton's parabola
}
\usage{
datagenParabola(n)
}
\arguments{
  \item{n}{
The sample size
}
}
\references{
Newton, M.A. (2009), 
Introducing the discussion paper by Sz?kely and Rizzo,
Annals of applied statistics, 3 (4), 1233-1235. 
}
\author{
Ruth Heller
}

\examples{
X = datagenParabola(50)
plot(X[1,],X[2,])
}