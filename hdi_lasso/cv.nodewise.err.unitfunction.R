cv.nodewise.err.unitfunction<-function (c, K, dataselects, x, lambdas, verbose, p) 
{
  if (verbose) {
    interesting.points <- round(c(1/4, 2/4, 3/4, 4/4) * p)
    names(interesting.points) <- c("25%", "50%", 
                                   "75%", "100%")
    if (c %in% interesting.points) {
      message("The expensive computation is now ", 
              names(interesting.points)[c == interesting.points], 
              " done")
    }
  }
  cv.nodewise.totalerr(c = c, K = K, dataselects = dataselects, 
                       x = x, lambdas = lambdas)
}