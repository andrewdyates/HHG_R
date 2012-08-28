\name{datagenW}
\alias{datagenW}
\title{
The W generating function
}
\description{
Newton's W
}
\usage{
datagenW(n)
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
X = datagenW(50)
plot(X[1,],X[2,])
}