
pvHHG<- function (Dx, Dy, M=1, alpha_hyp=NULL, alpha0=NULL, beta0=NULL, monte = 10000 , w_sum = 0, w_max = 2, score_id=0, epsilon=0.01) 
{
# score_id - decides which score will be used in the sequential testing (0 sum chi, 1 sum_like, 2 max chi, 3 max like)

#parameters of Sequential 
# all parameters till w_max (including) didnt change
#the next parameter (out2[[10+score_id]]) is the result of the Biometrika test on the real data (the true score).  
#M the overall number of hypotheses
# alpha_hyp is the alpha` you described (we are interested only in p_values below this), ie test H: p>alpha_hyp,  if unspecified then set it to 0.05/log(M) 
# epsilon - means we are trying to differentiate between alpha`*(1+epsilon) and alpha`(1-epsilon)
# alpha0 is the probability of not stopping early, even though p>alpha_hyp. 
# beta0 is the probability of stopping early and setting p=1, even though o<=alpha_hyp. This is the more severe error. 
#Return: 
# pv - is one if we hit B-threshold and otherwise is the p_val we got
# output_monte - how many iterations did we really run



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



# B_threshold is the threshold (log (beta0/(1-alpha0)))of the log likelihood ratio under which we stop. (natural logarithm), ie decide that the hypothesis is a null hypothesis
if (is.null(alpha_hyp))
alpha_hyp =0.05/ifelse(M==1,1,log(M))
if (is.null(alpha0))
alpha0 =0.05
if (is.null(beta0))
beta0 =min(0.01, 0.05/M)


B_threshold = log(beta0/(1-alpha0))
A_threshold = log((1-beta0)/alpha0)

output_pval=0
output_monte=0

out4=.C("Sequential",as.integer(n), as.integer(monte), as.double(as.matrix(Dx,nr=n,nc=n)),as.double(as.matrix(Dy,nr=n,nc=n)),as.integer(sum_chi_flag), as.integer(sum_like_flag),as.integer(max_chi_flag),
as.integer(max_like_flag),as.integer(w_sum),as.integer(w_max),as.integer(score_id),as.double(out2[[10+score_id]]), as.double(A_threshold), as.double(B_threshold),as.double(alpha_hyp),as.double(epsilon),as.double(output_pval),as.integer(output_monte))
pv = out4[[17]]
output_monte = out4[[18]]
return(list(pv=pv, output_monte = output_monte, A_threshold= A_threshold, B_threshold = B_threshold))
}

