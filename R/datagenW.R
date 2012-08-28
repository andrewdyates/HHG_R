datagenW <-
function(n){
x <- seq( -1, 1, length=n )

u <- x + runif(n)/3
v <-  4*( ( x^2 - 1/2 )^2 + runif(n)/500 )
return(rbind(u,v))
}
