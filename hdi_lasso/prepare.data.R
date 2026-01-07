prepare.data<-function (x, y, standardize, family) 
{
  x <- scale(x, center = TRUE, scale = standardize)
  dataset <- switch(family, gaussian = {
    list(x = x, y = y)
  }, {
    switch.family(x = x, y = y, family = family)
  })
  x <- scale(dataset$x, center = TRUE, scale = FALSE)
  y <- scale(dataset$y, scale = FALSE)
  y <- as.numeric(y)
  list(x = x, y = y)
}