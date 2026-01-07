Sum<-function (U, CovS) 
{
    if (all(abs(sum(U)) < 1e-20)) 
        pTsum <- 1
    else {
        a <- rep(1, length(U))
        Tsum <- sum(U)/(sqrt(as.numeric(t(a) %*% CovS %*% (a))))
        pTsum <- as.numeric(1 - pchisq(Tsum^2, 1))
    }
    pTsum
}