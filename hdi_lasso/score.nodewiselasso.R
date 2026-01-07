score.nodewiselasso<-function (x, wantTheta = FALSE, verbose = FALSE, lambdaseq = "quantile", 
          parallel = FALSE, ncores = 8, oldschool = FALSE, lambdatuningfactor = 1, 
          cv.verbose = FALSE, do.ZnZ = TRUE) 
{
  lambdas <- switch(lambdaseq, quantile = nodewise.getlambdasequence(x), 
                    linear = nodewise.getlambdasequence.old(x, verbose), 
                    stop("invalid 'lambdaseq': ", lambdaseq))
  if (verbose) {
    cat("Using the following lambda values:", lambdas, 
        "\n")
  }
  cvlambdas <- cv.nodewise.bestlambda(lambdas = lambdas, x = x, 
                                      parallel = parallel, ncores = ncores, oldschool = oldschool, 
                                      verbose = cv.verbose)
  if (verbose) {
    cat(paste("lambda.min is", cvlambdas$lambda.min), 
        "\n")
    cat(paste("lambda.1se is", cvlambdas$lambda.1se), 
        "\n")
  }
  if (do.ZnZ) {
    bestlambda <- improve.lambda.pick(x = x, parallel = parallel, 
                                      ncores = ncores, lambdas = lambdas, bestlambda = cvlambdas$lambda.min, 
                                      verbose = verbose)
    if (verbose) {
      cat("Doing Z&Z technique for picking lambda\n")
      cat("The new lambda is", bestlambda, "\n")
      cat("In comparison to the cross validation lambda, lambda = c * lambda_cv\n")
      cat("c=", bestlambda/cvlambdas$lambda.min, 
          "\n")
    }
  }
  else {
    if (lambdatuningfactor == "lambda.1se") {
      if (verbose) 
        cat("lambda.1se used for nodewise tuning\n")
      bestlambda <- cvlambdas$lambda.1se
    }
    else {
      if (verbose) 
        cat("lambdatuningfactor used is", lambdatuningfactor, 
            "\n")
      bestlambda <- cvlambdas$lambda.min * lambdatuningfactor
    }
  }
  if (verbose) {
    cat("Picked the best lambda:", bestlambda, "\n")
  }
  if (wantTheta) {
    out <- score.getThetaforlambda(x = x, lambda = bestlambda, 
                                   parallel = parallel, ncores = ncores, oldschool = TRUE, 
                                   verbose = verbose)
  }
  else {
    Z <- score.getZforlambda(x = x, lambda = bestlambda, 
                             parallel = parallel, ncores = ncores, oldschool = oldschool)
    out <- Z
  }
  return.out <- list(out = out, bestlambda = bestlambda)
  return(return.out)
}