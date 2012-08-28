\name{HHG2x2-package}
\alias{HHG2x2-package}
\alias{HHG2x2}
\docType{package}
\title{
A nonparametric test of independence between two random vectors based on distances.
}
\description{
A nonparametric test of independence between two random vectors, by using the pairwise distances between the obsrevations. 
Important: the computation of the test statistics is O(n^2log(n)), where n is the sample size. 
}
\details{
\tabular{ll}{
Package: \tab HHG\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2012-04-23\cr
License: \tab GPL-2
}

}
\author{
Ruth Heller <ruheller@post.tau.ac.il>,  Yair Heller <heller.yair@gmail.com>

Maintainer: Ruth Heller <ruheller@post.tau.ac.il>
}
\references{
Heller R., Heller Y., and Gorfine M. (2012). A consistent multivariate test of association based
on ranks of distances. \emph{arXiv:1201.3522v1}

}
\examples{

## a non-null univariate example

X = datagenCircle(50); 
plot(X[1,],X[2,])

Dx = as.matrix(dist((X[1,]),diag=TRUE,upper=TRUE))
Dy = as.matrix(dist((X[2,]),diag=TRUE,upper=TRUE))

n = dim(Dx)[1] ## the sample size

pvHHG(Dx,Dy,monte = 1000, M=1) 

## a null univariate example

X = datagen4indclouds(50); 

Dx = as.matrix(dist((X[1,]),diag=TRUE,upper=TRUE))
Dy = as.matrix(dist((X[2,]),diag=TRUE,upper=TRUE))

n = dim(Dx)[1] #the sample size

pvHHG(Dx,Dy,monte = 1000) 

## a multivariate example
library(MASS)

x = t(mvrnorm(50, rep(0,5),diag(1,5)))
y=log(x^2)
Dx = as.matrix(dist((t(x)),diag=TRUE,upper=TRUE))
Dy = as.matrix(dist((t(y)),diag=TRUE,upper=TRUE))
n = dim(Dx)[1] #the sample size

M=100
pvHHG(Dx,Dy,M, monte = 1000)
}
