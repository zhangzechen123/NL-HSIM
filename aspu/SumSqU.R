SumSqU<-function (U, CovS) 
{
    if (is.null(dim(CovS))) {
        Tscore <- sum(U^2/CovS)
        if (is.na(Tscore) || is.infinite(Tscore) || is.nan(Tscore)) 
            Tscore <- 0
        pTg1 <- as.numeric(1 - pchisq(Tscore, 1))
    }
    else {
        if (all(abs(U) < 1e-20)) 
            pTg1 <- 1
        else {
            Tg1 <- t(U) %*% U
            cr <- eigen(CovS, only.values = TRUE)$values
            alpha1 <- sum(cr * cr * cr)/sum(cr * cr)
            beta1 <- sum(cr) - (sum(cr * cr)^2)/(sum(cr * cr * 
                cr))
            d1 <- (sum(cr * cr)^3)/(sum(cr * cr * cr)^2)
            alpha1 <- as.double(alpha1)
            beta1 <- as.double(beta1)
            d1 <- as.double(d1)
            pTg1 <- as.numeric(1 - pchisq((Tg1 - beta1)/alpha1, 
                d1))
        }
    }
    return(pTg1)
}