# GPP响应、bootstrap和敏感性分析；加载时只定义函数。
# 敏感性分析使用drought_indices.R中的rolling_cor()。

# GPP分组描述统计。
group_summary <- function(d, col) {
  use <- is.finite(d$GPP) & !is.na(d[[col]])
  pieces <- split(d[use, ], d[[col]][use])
  do.call(rbind, lapply(names(pieces), function(nm) {
    x <- pieces[[nm]]$GPP
    data.frame(group = nm, n = length(x), mean_GPP = mean(x),
               median_GPP = median(x), sd_GPP = sd(x))
  }))
}

# 移动块重抽样与置信区间；依赖参数显式传入。
draw_blocks <- function(n, block_n) {
  starts <- sample.int(n - block_n + 1L, ceiling(n / block_n), replace = TRUE)
  unlist(lapply(starts, function(s) seq.int(s, s + block_n - 1L)))[seq_len(n)]
}

mean_difference <- function(x) {
  u <- x$GPP[x$coupling == "Synchronous" & !is.na(x$coupling)]
  v <- x$GPP[x$coupling == "Asynchronous" & !is.na(x$coupling)]
  if (!any(is.finite(u)) || !any(is.finite(v))) return(NA_real_)
  mean(u, na.rm = TRUE) - mean(v, na.rm = TRUE)
}

ci <- function(x, bootstrap_n) {
  x <- x[is.finite(x)]
  if (length(x) < bootstrap_n * 0.8) return(c(NA_real_, NA_real_))
  unname(quantile(x, c(0.025, 0.975)))
}

# GPP响应模型：控制温度/辐射异常、季节项和趋势。
# 没有两类耦合时段时只估计现有组；不通过改变阈值制造组间差异。
fit_gpp_response <- function(dat, bootstrap_n, block_n) {
  dat$coupling <- factor(dat$coupling, levels = c("Asynchronous", "Synchronous"))
  phase <- 2 * pi * as.integer(format(dat$date, "%j")) / 365.25
  dat$season_sin <- sin(phase)
  dat$season_cos <- cos(phase)
  dat$time_year <- as.numeric(dat$date - min(dat$date)) / 365.25
  formula_gpp <- GPP ~ coupling * (D_VPD + D_SM1) + Z_Tavg + Z_Rs +
    season_sin + season_cos + time_year
  observed_groups <- unique(as.character(dat$coupling[is.finite(dat$GPP) & !is.na(dat$coupling)]))
  two_groups <- length(observed_groups) == 2L
  if (!two_groups) {
    # 没有两类样本时不强行改变阈值，改为估计现有组的响应。
    formula_gpp <- GPP ~ D_VPD + D_SM1 + Z_Tavg + Z_Rs + season_sin + season_cos + time_year
  }
  model_cols <- all.vars(formula_gpp)
  keep <- complete.cases(dat[model_cols])
  fit <- lm(formula_gpp, data = dat[keep, ])
  if (anyNA(coef(fit))) stop("回归存在不可估计系数，请检查分组样本量")


  # 从首个到末个有效GPP之间的完整时间轴重抽样，缺测不会被压缩连接。
  span <- seq.int(min(which(is.finite(dat$GPP))), max(which(is.finite(dat$GPP))))
  boot_data <- dat[span, ]
  if (any(as.numeric(diff(boot_data$date)) > 8)) stop("GPP研究期存在时间缺口")
  n <- nrow(boot_data)
  block_n <- min(block_n, n)


  boot_coef <- matrix(NA_real_, bootstrap_n, length(coef(fit)),
                      dimnames = list(NULL, names(coef(fit))))
  boot_delta <- rep(NA_real_, bootstrap_n)
  for (i in seq_len(bootstrap_n)) {
    d <- boot_data[draw_blocks(n, block_n), ]
    boot_delta[i] <- mean_difference(d)
    cb <- tryCatch(coef(lm(formula_gpp, data = d)), error = function(e) NULL)
    if (!is.null(cb)) boot_coef[i, intersect(names(cb), colnames(boot_coef))] <-
      cb[intersect(names(cb), colnames(boot_coef))]
  }

  intervals <- t(apply(boot_coef, 2, ci, bootstrap_n = bootstrap_n))
  coefficients <- data.frame(term = names(coef(fit)), estimate = unname(coef(fit)),
                              lower_95 = intervals[, 1], upper_95 = intervals[, 2],
                              bootstrap_valid = colSums(is.finite(boot_coef)))

  delta_ci <- ci(boot_delta, bootstrap_n)
  difference <- data.frame(comparison = "Synchronous minus Asynchronous",
                        GPP_mean_difference = mean_difference(boot_data),
                        lower_95 = delta_ci[1], upper_95 = delta_ci[2])

  # 将交互项转换为每一组的干旱响应斜率，便于解释。
  slopes <- do.call(rbind, lapply(c("D_VPD", "D_SM1"), function(v) {
    interaction <- paste0("couplingSynchronous:", v)
    do.call(rbind, lapply(c("Asynchronous", "Synchronous"), function(g) {
      if (!g %in% observed_groups) {
        return(data.frame(index = v, coupling = g, slope = NA_real_,
                           lower_95 = NA_real_, upper_95 = NA_real_))
      }
      value <- coef(fit)[v]
      samples <- boot_coef[, v]
      if (two_groups && g == "Synchronous") {
        value <- value + coef(fit)[interaction]
        samples <- samples + boot_coef[, interaction]
      }
      limits <- ci(samples, bootstrap_n)
      data.frame(index = v, coupling = g, slope = unname(value),
                 lower_95 = limits[1], upper_95 = limits[2])
    }))
  }))
  list(data = dat, model = fit, coefficients = coefficients,
       mean_difference = difference, slopes = slopes, two_groups = two_groups,
       block_n = block_n)
}

# 预设窗口和各土层的敏感性分析。
window_layer_sensitivity <- function(dat, cor_cut) {
  sensitivity <- do.call(rbind, lapply(c(15L, 23L, 46L), function(w) {
    do.call(rbind, lapply(paste0("SM", 1:4), function(v) {
      r <- rolling_cor(dat$D_VPD, dat[[paste0("D_", v)]], dat$date, w)
      d <- dat
      d$coupling <- ifelse(is.na(r), NA_character_,
                           ifelse(r >= cor_cut, "Synchronous", "Asynchronous"))
      z <- group_summary(d, "coupling")
      z$window_periods <- w
      z$soil_layer <- v
      z
    }))
  }))
  sensitivity
}

