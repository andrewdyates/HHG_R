\name{pvHHG.slow}
\alias{pvHHG.slow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the new tests p-values when the 2 by 2 tables are summarized either by Pearson's statistics or the likelihood ratio statistic . 
}
\description{
The p-value is based on the null distribution approximated by sampling random permutations of the labels of the y vector.
Important: the computation of the test statistics is O(n^3), where n is the sample size. Therefore, it may take a while. For example, with 100000 monte-carlo samples and n=50 the test takes 3 minutes to run on my computer. 
}
\usage{
pvHHG.slow(Dx,Dy,monte=10000, w_sum=0, w_max = 2)
}
\arguments{
  \item{Dx}{
The pairwise distances in X,  .
}

  \item{Dy}{
The pairwise distances in Y.
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
}
\details{
Compute the 4 new test p-values. 
}
\value{
  \item{pv.sum.chisq }{The pvalue of the test statitisic in Biometrika: sum of Pearson test statistics across all 2x2 tables}
 \item{pv.sum.LR }{The pvalue of the alternative discussed in Biometrika: sum of LR test statistics across all 2x2 tables}
  \item{pv.max.chisq }{The pvalue of the test statitisic that is the maximum of the Pearson test statistics across all 2x2 tables}
 \item{pv.max.LR }{The pvalue of the statistic that is the maximum  of LR test statistic across all 2x2 tables}

}
\references{
Heller R., Heller Y., and Gorfine M., "A consistent multivariate test of association based
on ranks of distances", arXiv:1201.3522v1 
}
\author{
Ruth Heller and Yair Heller
}

