get.clusterGroupTest.function<-function (group.testing.function, x) 
{
  stopifnot(is.function(group.testing.function))
  function(hcloutput, dist = as.dist(1 - abs(cor(x))), alpha = 0.05, 
           method = "average", conservative = TRUE) {
    hh <- if (missing(hcloutput)) 
      hclust(dist, method = method)
    else hcloutput
    clusterextractfunction <- function(nclust) {
      cutree(hh, k = nclust)
    }
    pvalfunction <- function(clusters.to.test, nclust) {
      mapply(group.testing.function, conservative = conservative, 
             group = split(clusters.to.test, col(clusters.to.test)))
    }
    structure(c(calculate.pvalue.for.cluster(hh = hh, p = ncol(x), 
                                             pvalfunction = pvalfunction, alpha = alpha, clusterextractfunction = clusterextractfunction), 
                method = "clusterGroupTest"), class = c("clusterGroupTest", 
                                                        "hdi"))
  }
}