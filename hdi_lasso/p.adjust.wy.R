p.adjust.wy<-function (cov, pval, N = 10000) 
{ # 检查 cov 和 pval 是否为标量
  if (is.null(dim(cov)) || (length(cov) == 1 && length(pval) == 1)) {
    # 对于标量情况，直接返回 pval
    return(pval)
  }
  zz <- mvrnorm(N, rep(0, ncol(cov)), cov)
  zz2 <- scale(zz, center = FALSE, scale = sqrt(diag(cov)))
  Gz <- apply(2 * pnorm(abs(zz2), lower.tail = FALSE), 1, min)
  # 使用经验分布函数计算 p 值
  return(ecdf(Gz)(pval))
}