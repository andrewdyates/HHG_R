\name{pvHHG}
\alias{pvHHG}
\title{
A nonparametric test of independence between two random vectors based on distances. 
}
\description{
A nonparametric test of independence between two random vectors, by using the pairwise distances between the obsrevations. 
Important: the computation of the test statistics is O(n^2log(n)), where n is the sample size. Therefore, it may take a while. 
}
\usage{
pvHHG(Dx, Dy, M=1, alpha_hyp=NULL, alpha0=NULL, beta0=NULL, monte = 10000 , w_sum = 0, w_max = 2, score_id=0, epsilon=0.01) 

}

\arguments{
  \item{Dx}{
The pairwise distances in X .
}

  \item{Dy}{
The pairwise distances in Y.
}

  \item{M}{
The total number of hypotheses in the family of hypotheses simultaneously tested.
}




  \item{alpha_hyp}{
 the value above which we can safely set the estimated  pvalue as 1 since this p-value will not lead to discoveries when faced with multiplicity of hypotheses;
 if unspecified, the vlaue is set to 0.05/log(M) if M>1, 0.05 for M=1. 
}

\item{alpha0}{
The probability of deciding that p<= alpha_hyp even though p>alpha_hyp; if unspecified, the value is set to 0.05. 
}


\item{beta0}{
The probability of deciding that p> alpha_hyp even though p<=alpha_hyp; if unspecified, the value is set to 0.05/M. 
}

  \item{monte}{
The number of random permutations of the labels of the Y vector used for approximating the null distribution of the test statistic.
}
  \item{w_sum}{
	minimum number of elements expected in cell for chi-square score  of the 2x2 table to be added to the score that is the sum of 2x2 table chi-squared statistics.
}
  \item{w_max}{
	minimum number of elements expected in cell for chi-sqaure score of the 2x2 table to be considered for the score that is the maximum of 2x2 table chi-squared statistics. 

}
  \item{score_id}{
	indicator of which score will be used in the sequential testing (0 = chi-square score, 1 = likelihood ratio score, 2 = chi-squared-max score, 3 = likelihood ratio-max score)
}
  \item{epsilon}{
in the sequential testing procedure, it means we are trying to differentiate between p_val_threshold*(1+epsilon) and p_val_threshold*(1-epsilon)
}
}

\details{ A nonparametric test of independence between two random vectors, by using the pairwise distances between the obsrevations. For details see Heller et al. (2012).
When simultaneous inference on M such tests is desired, and only a handful R of these tests is expected to lead to discoveries, we recommend calling the test with parameter p_val_threshold = R*0.05/M for each test 
for the computational purposes. The computation of the test statistics is O(n^2log(n)), where n is the sample size. 
We implement Wald's sequential test, as suggested by Fay et al. (2007), to reduce the cost of O(monte*n^2log(n))
to a smaller number. When faced with multiplicity, then necessary number of iterations monte may be quite large. For example, if it is desired to control the FWER at level 0.05 with a Bonferroni correction for multiplicity, 
the minimum number of iterations to establish significance of a p-value is monte = 20*M. However, if 
the hypothesis is deemed uninteresting after a smaller number of iteration, the function may be exited and the p-value estimated to be 1. Therefore, in practice when only a handful of hypotheses are expected to be 
non-null the number of iterations in total, for all hypotheses, may be actually lmuch ower than 20*M*M. The threshold for decideding which hypothesis is interesting depends on the 
total number of hypotheses simultaneously examined. For a Bonferroni correction, this threshold is alpha/M, where alpha is the desired FWER level (say alpha=0.05). 
For the less conservative procedure of Benjamini & Hochberg (1995), this threshold is M1*q/M, where q is the desired FDR level (say q=0.05), and M1 is the unknwon number of true non-null hypothesese (ie dependent associations). 
Although M1 is unknown, it may conservatively be estimated by the investigator (e.g. if at most 0.02 of the hypotheses are expected to be non-null, set M1 = 0.02xM).  
}
\value{
  \item{pv }{The pvalue of the test; it is one if we hit the  B-threshold and otherwise it is the p-value estimated from "monte" number of iterations}
 \item{monte }{The actual number of monte-carlo iterations used}

}
\references{
Heller R., Heller Y., and Gorfine M. (2012). A consistent multivariate test of association based
on ranks of distances. \emph{arXiv:1201.3522v1}

Fay, M., and Kim., H., and Hachey, M. (2007) . On Using Truncated Sequential Probability Ratio Test Boundaries for Monte Carlo Implementation of Hypothesis Tests. \emph{J Comput Graph Stat}, 16(4): 946-967.

Benjamini, Y., and Hochberg, Y. (1995).Controlling the false discovery rate: a practical and powerful approach to multiple testing.
  \emph{Journal of the Royal Statistical Society Series} B, \bold{57}, 289--300.
  
}
\author{
Ruth Heller and Yair Heller
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
