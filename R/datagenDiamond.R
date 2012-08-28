datagenDiamond <-
function(n){
x <- runif(n, min=(-1), max=1 )
y <- runif(n, min=(-1), max=1 )

theta <- -pi/4
rr <- rbind( c(cos(theta), -sin(theta) ),
             c( sin(theta), cos(theta) ) )
tmp <- cbind( x, y ) %*% rr
u <- tmp[,1]
v <-  tmp[,2]
return(rbind(u,v))
}
