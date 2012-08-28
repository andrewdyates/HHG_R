Modified HHG implementation to directly return the T statistics.

Compile:

    R CMD INSTALL HHG2x2 -l /path/to/HHG_R/

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
