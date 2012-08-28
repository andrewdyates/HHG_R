datagen2Parabolas <-
function(n){
x <- seq(-1,1, length=n )
y <- (x ^2 + runif(n)/2 )*( sample( c(-1,1), size=n, replace=T ) )
return(rbind(x,y))
}
