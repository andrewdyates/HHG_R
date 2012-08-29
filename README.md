Modified HHG implementation; returns the T statistic rather than merely the p-value.

Maximum values (log base _e_):

  * sum_chisquared: `(n)*(n-1)*(n-2)`
  * sum_lr: `(n)*(n-1)*(n-2) / log(2)`
  * max_chisquared: `n-2`
  * max_chisquared: `n-2 / log(2)`

Compile:

    cd /path/to/HHG_R
    cd ..
    R CMD INSTALL HHG_R

Run in R:

    > library("HHG2x2", lib.loc="/path/to/HHG_R")
    > X = datagenCircle(50);
    > Dx = as.matrix(dist((X[1,]),diag=TRUE,upper=TRUE))
    > Dy = as.matrix(dist((X[2,]),diag=TRUE,upper=TRUE))
    > myHHG(Dx,Dy)
    > pvHHG(Dx,Dy)

Output of myHHG:

    $sum_chisquared
    [1] 5081.871
    
    $sum_lr
    [1] 2815.306
    
    $max_chisquared
    [1] 16.45584
    
    $max_lr
    [1] 10.03125


See test_hhgR.py for python integration examples.

Errors:
    rpy2.rinterface.RRuntimeError: Error: package ‘HHG2x2’ was built for x86_64-apple-darwin10.8.0

    
Delete HHG2x2.o and HHG2x2.so and rebuild.


