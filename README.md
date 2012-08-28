Compile:

    R CMD INSTALL HHG2x2 -l /Users/z/Desktop/HHG2x2/

Run in R:

    > library("HHG2x2", lib.loc="/Users/z/Desktop/HHG2x2")
    > X = datagen4indclouds(50); 
    > Dx = as.matrix(dist((X[1,]),diag=TRUE,upper=TRUE))
    > Dy = as.matrix(dist((X[2,]),diag=TRUE,upper=TRUE))
    > myHHG(Dx,Dy)

Output:

    $sum_chisquared
    [1] 5081.871
    
    $sum_lr
    [1] 2815.306
    
    $max_chisquared
    [1] 16.45584
    
    $max_lr
    [1] 10.03125
