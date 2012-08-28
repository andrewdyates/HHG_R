datagen4indclouds <-
function(n){
dx <- rnorm(n)/3
dy <- rnorm(n)/3
cx <- sample( c(-1,1), size=n, replace=T )
cy <- sample( c(-1,1), size=n, replace=T )
u <- cx + dx
v <- cy + dy 
return(rbind(u,v))
}
