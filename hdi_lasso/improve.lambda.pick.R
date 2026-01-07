improve.lambda.pick<-function (x, parallel, ncores, lambdas, bestlambda, verbose) 
{
  lambdas <- sort(lambdas, decreasing = TRUE)
  M <- calcM(x = x, lambdas = lambdas, parallel = parallel, 
             ncores = ncores)
  Mcv <- M[which(lambdas == bestlambda)]
  if (length(which(M < 1.25 * Mcv)) > 0) {
    lambdapick <- min(lambdas[which(M < 1.25 * Mcv)])
  }
  else {
    if (verbose) 
      cat("no better lambdapick found\n")
    lambdapick <- bestlambda
  }
  if (max(which(M < 1.25 * Mcv)) < length(lambdas)) {
    if (verbose) 
      cat("doing a second step of discretisation of the lambda space to improve the lambda pick\n")
    lambda.number <- max(which(M < 1.25 * Mcv))
    newlambdas <- seq(lambdas[lambda.number], lambdas[lambda.number + 
                                                        1], (lambdas[lambda.number + 1] - lambdas[lambda.number])/100)
    newlambdas <- sort(newlambdas, decreasing = TRUE)
    M2 <- calcM(x = x, lambdas = newlambdas, parallel = parallel, 
                ncores = ncores)
    if (length(which(M2 < 1.25 * Mcv)) > 0) {
      evenbetterlambdapick <- min(newlambdas[which(M2 < 
                                                     1.25 * Mcv)])
    }
    else {
      if (verbose) 
        cat("no -even- better lambdapick found\n")
      evenbetterlambdapick <- lambdapick
    }
    if (is.infinite(evenbetterlambdapick)) {
      if (verbose) {
        cat("hmmmm the better lambda pick after the second step of discretisation is Inf\n")
        cat("M2 is\n")
        cat(M2, "\n")
        cat("Mcv is\n")
        cat(Mcv, "\n")
        cat("and which(M2 < 1.25* Mcv) is\n")
        cat(which(M2 < 1.25 * Mcv), "\n")
      }
    }
    else {
      lambdapick <- evenbetterlambdapick
    }
  }
  return(lambdapick)
}