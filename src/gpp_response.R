# GPP分组均值与Pearson相关分析；加载时只定义函数。

paired_correlation <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  r <- if (n >= 3L && sd(x[ok]) > 0 && sd(y[ok]) > 0) cor(x[ok], y[ok]) else NA_real_
  list(n = n, r = r)
}

group_summary <- function(dat, col, groups) {
  do.call(rbind, lapply(groups, function(group) {
    d <- dat[!is.na(dat[[col]]) & dat[[col]] == group, ]
    gpp <- d$GPP[is.finite(d$GPP)]
    vpd <- paired_correlation(d$GPP, d$D_VPD)
    sm <- paired_correlation(d$GPP, d$D_SM1)
    data.frame(grouping = col, group = group, n_GPP = length(gpp),
               mean_GPP = if (length(gpp)) mean(gpp) else NA_real_,
               n_VPD = vpd$n, r_VPD = vpd$r, n_SM1 = sm$n, r_SM1 = sm$r)
  }))
}

analyze_gpp_response <- function(dat) {
  rbind(group_summary(dat, "coupling", c("Asynchronous", "Synchronous")),
        group_summary(dat, "drought_state", c("Both", "VPD_only", "SM_only", "Neither")))
}
