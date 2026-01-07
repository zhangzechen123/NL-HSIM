do.initial.fit<-function (x, y, initial.lasso.method = c("scaled lasso", 
                                           "cv lasso"), lambda, verbose = FALSE) 
{
  no.lambda.given <- missing(lambda) || is.null(lambda)
  if (!no.lambda.given && verbose) {
    cat("A value for lambda was provided to the do.initial.fit function,\n")
    cat("the initial.lasso.method option was therefore ignored.\n")
    cat("We now do a lasso with the tuning parameter instead of a self-tuning procedure.\n")
  }
  if (no.lambda.given) {
    switch(initial.lasso.method, `scaled lasso` = {
      scaledlassofit <- scalreg(X = x, y = y)
      lambda <- NULL
    }, `cv lasso` = {
      glmnetfit <- cv.glmnet(x = x, y = y)
      lambda <- glmnetfit$lambda.1se
    }, {
      stop("Not sure what lasso.method you want me to use for the initial fit. The only options for the moment are: 1)scaled lasso 2)cvlasso")
    })
  }
  else {
    glmnetfit <- glmnet(x = x, y = y)
  }
  if (no.lambda.given && identical(initial.lasso.method, "scaled lasso")) {
    intercept <- 0
    betalasso <- scaledlassofit$coefficients
    sigmahat <- scaledlassofit$hsigma
    residual.vector <- y - x %*% betalasso
  }
  else {
    if ((nrow(x) - sum(as.vector(coef(glmnetfit, s = lambda)) != 
                       0)) <= 0) {
      message("Refitting using cross validation: your lambda used all degrees of freedom")
      glmnetfit <- cv.glmnet(x = x, y = y)
      lambda <- glmnetfit$lambda.1se
    }
    intercept <- coef(glmnetfit, s = lambda)[1]
    betalasso <- as.vector(coef(glmnetfit, s = lambda))[-1]
    residual.vector <- y - predict(glmnetfit, newx = x, s = lambda)
    sigmahat <- sqrt(sum((residual.vector)^2)/(nrow(x) - 
                                                 sum(as.vector(coef(glmnetfit, s = lambda)) != 0)))
  }
  list(betalasso = betalasso, sigmahat = sigmahat, intercept = intercept, 
       lambda = lambda)
}