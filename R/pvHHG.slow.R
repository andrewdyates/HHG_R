pvHHG.slow <-
function (Dx, Dy, monte = 10000, w_sum = 0, w_max = 2) 
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
    v_sc = seq(0, length = monte)
    v_sl = seq(0, length = monte)
    v_mc = seq(0, length = monte)
    v_ml = seq(0, length = monte)
    cumulative = seq(0, length = monte)
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
    out3 = .C("Create_Biometrika_Distribution", as.integer(n), 
        as.integer(monte), as.double(as.matrix(Dx, nr = n, nc = n)), 
        as.double(as.matrix(Dy, nr = n, nc = n)), as.integer(sum_chi_flag), 
        as.integer(sum_like_flag), as.integer(max_chi_flag), 
        as.integer(max_like_flag), as.integer(w_sum), as.integer(w_max), 
        as.double(as.vector(v_sc)), as.double(as.vector(v_sl)), 
        as.double(as.vector(v_mc)), as.double(as.vector(v_ml)), 
        as.double(as.vector(cumulative)))
    out.sum.chisquared.dist = out3[[11]]
    out.sum.lr.dist = out3[[12]]
    out.max.chisquared.dist = out3[[13]]
    out.max.lr.dist = out3[[14]]
    out.pdist = out3[[15]]
    pv.sum.chisq = 1 - out.pdist[sum(out.sum.chisquared.dist <= 
        out.sum.chisquared.stat)]
    pv.sum.LR = 1 - out.pdist[sum(out.sum.lr.dist <= out.sum.lr.stat)]
    pv.max.chisq = 1 - out.pdist[sum(out.max.chisquared.dist <= 
        out.max.chisquared.stat)]
    pv.max.LR = 1 - out.pdist[sum(out.max.lr.dist <= out.max.lr.stat)]
    return(list(pv.sum.chisq = pv.sum.chisq, pv.sum.LR = pv.sum.LR, 
        pv.max.chisq = pv.max.chisq, pv.max.LR = pv.max.LR))
}
