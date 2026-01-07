PowerUniv<-function (U, V) 
{
    n <- dim(V)[1]
    x <- as.numeric(max(abs(U)))
    TER <- as.numeric(1 - pmvnorm(lower = c(rep(-x, n)), upper = c(rep(x, 
        n)), mean = c(rep(0, n)), sigma = V))
    return(TER)
}