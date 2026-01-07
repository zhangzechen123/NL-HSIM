UminPd<-function (U, CovS) 
{
    if (is.null(dim(CovS))) {
        Tu <- sum(U^2/CovS)
        if (is.na(Tu) || is.infinite(Tu) || is.nan(Tu)) 
            Tu <- 0
        pTu <- as.numeric(1 - pchisq(Tu, 1))
    }
    else {
        Tu <- as.vector(abs(U)/(sqrt(diag(CovS)) + 1e-20))
        k <- length(U)
        V <- matrix(0, nrow = k, ncol = k)
        for (i in 1:k) {
            for (j in 1:k) {
                if (abs(CovS[i, j]) > 1e-20) 
                  V[i, j] <- CovS[i, j]/sqrt(CovS[i, i] * CovS[j, 
                    j])
                else V[i, j] <- 1e-20
            }
        }
        pTu <- as.numeric(PowerUniv(Tu, V))
    }
    pTu
}