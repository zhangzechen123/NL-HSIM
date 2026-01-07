ART.A <- function(P, k, L) { # "ART-A"
    wgt <- rep(1,k)
    z <- P
    z[1] <- ( 1 - P[1] )^L
    for(j in 2:L) z[j] <- ((1-P[j]) / (1-P[j-1]))^((L-(j-1)))
    p <- (1-z)[1:k]
    k = length(p)
    sumZ <- rep(0, k)
    y <- qnorm(p)
    z <- y
    gz <- z[1] * wgt[1]
    sumZ[1] <- gz
    for(i in 2:k) {
        gz <- p[ 1 : i ]
        for(j in 1:i) gz[j] <- z[j] * wgt[j]
        sumZ[i] <- sum(gz)
    }
    Lo = diag(k); Lo[lower.tri(Lo)] <- 1
    pSg <- Lo %*% diag(wgt[1:k]^2) %*% t(Lo)
    pCr <<- cov2cor(pSg)
    sZ <- sumZ
    for(i in 1:k) {
        sZ[i] <- sumZ[i] / sqrt(diag(pSg))[i]
    }
    sZ<<-sZ
    ppZ <- pmvnorm(lower = rep(-Inf,k), upper = rep(max(sZ), k), sigma = pCr)[1]
    c(ppZ, which.min(sZ))
}
