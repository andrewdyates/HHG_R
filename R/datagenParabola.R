datagenParabola <-
function(n){
x <- seq(-1,1, length=n )
y <- (x ^2 + runif(n))/2
return(rbind(x,y))
}
