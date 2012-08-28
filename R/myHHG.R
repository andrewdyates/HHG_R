


myHHG<- function (Dx, Dy, w_sum = 0, w_max = 2) 
{
    n = dim(Dx)[1]
    cutoff = 0
    sum_chi_flag = 1
    sum_like_flag = 1
    max_chi_flag = 1
    max_like_flag = 1
    sc = 0
    sl = 0
    mc = 0
    ml = 0
    out2 = .C("Biometrika", as.integer(n), as.double(as.matrix(Dx, 
        nr = n, nc = n)), as.double(as.matrix(Dy, nr = n, nc = n)), 
        as.integer(sum_chi_flag), as.integer(sum_like_flag), 
        as.integer(max_chi_flag), as.integer(max_like_flag), 
        as.integer(w_sum), as.integer(w_max), as.double(sc), 
        as.double(sl), as.double(mc), as.double(ml))
    out.sum.chisquared.stat = out2[[10]]
    out.sum.lr.stat = out2[[11]]
    out.max.chisquared.stat = out2[[12]]
    out.max.lr.stat = out2[[13]]
    normed_sum_chisquared = out2[[10]]/(n*(n-1)*(n-2))

    return(list(sum_chisquared = out2[[10]], sum_lr=out2[[11]], max_chisquared=out2[[12]], max_lr=out2[[13]], normed_sum_chisquared = out2[[10]]/n/(n-1)/(n-2) ))

}


