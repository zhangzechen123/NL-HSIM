initial.estimator<-function (betainit, x, y, sigma) 
{
  if (!((is.numeric(betainit) && length(betainit) == ncol(x)) || 
        (betainit %in% c("scaled lasso", "cv lasso")))) {
    stop("The betainit argument needs to be either a vector of length ncol(x) or one of 'scaled lasso' or 'cv lasso'")
  }
  warning.sigma.message <- function() {
    warning("Overriding the error variance estimate with your own value. The initial estimate implies an error variance estimate and if they don't correspond the testing might not be correct anymore.")
  }
  lambda <- NULL
  if (is.numeric(betainit)) {
    beta.lasso <- betainit
    if (is.null(sigma)) {
      stop("Not sure what variance estimate to use here")
      residual.vector <- y - x %*% betainit
      sigmahat <- sqrt(sum((residual.vector)^2)/(nrow(x) - 
                                                   sum(betainit != 0)))
    }
    else {
      warning.sigma.message()
      sigmahat <- sigma
    }
  }
  else {
    initial.fit <- do.initial.fit(x = x, y = y, initial.lasso.method = betainit)
    beta.lasso <- initial.fit$betalasso
    lambda <- initial.fit$lambda
    if (is.null(sigma)) {
      sigmahat <- initial.fit$sigmahat
    }
    else {
      warning.sigma.message()
      sigmahat <- sigma
    }
  }
  list(beta.lasso = beta.lasso, sigmahat = sigmahat, lambda = lambda)
}