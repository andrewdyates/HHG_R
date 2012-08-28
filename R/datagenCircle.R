datagenCircle <-
function(n){
x <- seq( -1, 1, length=n )
u <- sin( x*pi ) + rnorm( n )/8
v <- cos( x*pi ) + rnorm( n )/8
return(rbind(u,v))
}
