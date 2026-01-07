
aSPUd2<-function (Y, X, cov = NULL, model = c("gaussian", "binomial")) 
{
  model = match.arg(model)
  n <- length(Y)
  if (is.null(X) && length(X) > 0) 
    X = as.matrix(X, ncol = 1)
  k <- ncol(X)
  if (is.null(cov)) {
    XUs <- Xg <- X
    r <- Y - mean(Y)
    U <- as.vector(t(Xg) %*% r)
    Xbar <- apply(Xg, 2, mean)
    Xgb <- Xg
    for (i in 1:nrow(Xg)) Xgb[i, ] <- Xg[i, ] - Xbar
    if (model == "binomial") {
      CovS <- mean(Y) * (1 - mean(Y)) * (t(Xgb) %*% Xgb)
    }
    else {
      CovS <- var(Y) * (t(Xgb) %*% Xgb)
    }
  }
  else {
    tdat1 <- data.frame(trait = Y, cov)
    fit1 <- glm(trait ~ ., family = model, data = tdat1)
    pis <- fitted.values(fit1)
    XUs <- matrix(0, nrow = n, ncol = k)
    Xmus = X
    for (i in 1:k) {
      tdat2 <- data.frame(X1 = X[, i], cov)
      fit2 <- glm(X1 ~ ., data = tdat2)
      Xmus[, i] <- fitted.values(fit2)
      XUs[, i] <- (X[, i] - Xmus[, i])
    }
    r <- Y - pis
    U <- t(XUs) %*% r
    Xg <- X
    Xbar <- apply(Xg, 2, mean)
    Xgb <- Xg
    for (i in 1:nrow(Xg)) Xgb[i, ] <- Xg[i, ] - Xbar
    if (model == "binomial") {
      CovS <- mean(pis * (1 - pis)) * (t(Xgb) %*% Xgb)
    }
    else {
      CovS <- var(r) * (t(Xgb) %*% Xgb)
    }
  }
  p1 = SumSqU(U, CovS)
  p2 = Sum(U, CovS)
  p3 = UminPd(U, CovS)
  mp = min(p1, p2, p3)
  p = 1 - (1 - mp)^3
  pvs = c(p1, p2, p3, p)
  names(pvs) = c("Psum", "Pssu", "PminP", "PaSPU")
  pvs
}