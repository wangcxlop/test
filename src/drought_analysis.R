# 干旱分析函数库：只定义函数，加载时不读取数据、不输出文件。
# 仅使用base R；GPP已为异常值，不再重复去季节化。

# 文件读写与数据质量检查 ------------------------------------------------
save_table <- function(x, name, out_dir) {
  write.csv(x, file.path(out_dir, name), row.names = FALSE,
            fileEncoding = "UTF-8", na = "")
}

read_input <- function(name, required, data_dir) {
  path <- file.path(data_dir, name)
  if (!file.exists(path)) stop("找不到文件：", path)
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(required %in% names(x))) stop(name, "缺少必要列")
  x$date <- as.Date(x$date)
  if (anyNA(x$date) || anyNA(x$site)) stop(name, "存在无效日期或站点")
  if (anyDuplicated(x[c("site", "date")])) stop(name, "存在重复站点日期")
  x[order(x$site, x$date), ]
}

read_drought_inputs <- function(data_dir) {
  met <- read_input(data_dir = data_dir, "ERA5L_TWSA.csv",
                    c("site", "date", "VPD", paste0("SM", 1:4), "Tavg", "Rs"))
  gpp <- read_input(data_dir = data_dir, "GPP_anomaly.csv", c("site", "date", "GPP"))
  if (length(unique(met$site)) != 1L || length(unique(gpp$site)) != 1L ||
      unique(met$site) != unique(gpp$site)) stop("本脚本要求两份文件为同一个站点")
  vars <- c("VPD", paste0("SM", 1:4), "Tavg", "Rs")
  if (!all(vapply(met[vars], is.numeric, logical(1))) || !is.numeric(gpp$GPP)) {
    stop("分析变量必须为数值列")
  }
  if (any(met$VPD < 0, na.rm = TRUE) ||
      any(as.matrix(met[paste0("SM", 1:4)]) < 0, na.rm = TRUE) ||
      any(as.matrix(met[paste0("SM", 1:4)]) > 1, na.rm = TRUE)) {
    stop("VPD或土壤体积含水量超出合理范围，请核对单位")
  }
  quality <- do.call(rbind, lapply(list(met = met, gpp = gpp), function(x) {
    data.frame(variable = names(x), missing = vapply(x, function(z) sum(is.na(z)), 0L))
  }))
  quality$dataset <- rep(c("met", "gpp"), c(ncol(met), ncol(gpp)))
  list(met = met, gpp = gpp, quality = quality, vars = vars)
}

# 按年内46个8天时段，以全期数据计算季节性z分数（不是SPI/SPEI）。
# SM指数取负号，统一为数值越大越干；全期基准用于回顾分析。
seasonal_z <- function(x, slot) {
  z <- rep(NA_real_, length(x))
  for (s in unique(slot)) {
    ids <- which(slot == s & is.finite(x))
    if (length(ids) >= 5L && sd(x[ids]) > 0) {
      z[ids] <- (x[ids] - mean(x[ids])) / sd(x[ids])
    }
  }
  z
}

# 连续干旱事件；强度为时间加权均值，烈度为指数乘天数的累计值。
event_table <- function(dates, index, days, label, threshold) {
  dry <- is.finite(index) & index >= threshold
  gap <- c(TRUE, as.numeric(diff(dates)) > 8)
  starts <- which(dry & (gap | !c(FALSE, head(dry, -1))))
  empty <- data.frame(index = character(), start = as.Date(character()),
                      end = as.Date(character()), n_periods = integer(),
                      duration_days = numeric(), intensity = numeric(),
                      peak = numeric(), severity_z_days = numeric(),
                      excess_severity_z_days = numeric(), boundary_or_gap = logical())
  if (!length(starts)) return(empty)
  do.call(rbind, lapply(starts, function(i) {
    j <- i
    while (j < length(dry) && dry[j + 1L] && !gap[j + 1L]) j <- j + 1L
    ids <- i:j
    censored <- i == 1L || j == length(dry)
    if (i > 1L) censored <- censored || gap[i] || !is.finite(index[i - 1L])
    if (j < length(dry)) censored <- censored || gap[j + 1L] || !is.finite(index[j + 1L])
    data.frame(index = label, start = dates[i], end = dates[j] + days[j] - 1,
               n_periods = length(ids), duration_days = sum(days[ids]),
               intensity = weighted.mean(index[ids], days[ids]), peak = max(index[ids]),
               severity_z_days = sum(index[ids] * days[ids]),
               excess_severity_z_days = sum((index[ids] - threshold) * days[ids]),
               boundary_or_gap = censored)
  }))
}

# 向后滑动相关；有效配对不少于80%，不跨越超过8天的日期缺口。
rolling_cor <- function(x, y, dates, width) {
  ans <- rep(NA_real_, length(x))
  if (length(x) < width) return(ans)
  for (i in seq.int(width, length(x))) {
    ids <- seq.int(i - width + 1L, i)
    if (any(as.numeric(diff(dates[ids])) > 8)) next
    ok <- is.finite(x[ids]) & is.finite(y[ids])
    if (sum(ok) >= ceiling(width * 0.8) &&
        sd(x[ids][ok]) > 0 && sd(y[ids][ok]) > 0) {
      ans[i] <- cor(x[ids][ok], y[ids][ok])
    }
  }
  ans
}

prepare_drought_indices <- function(met, vars, threshold, window_n, cor_cut) {
  met$slot <- (as.integer(format(met$date, "%j")) - 1L) %/% 8L + 1L

  for (v in vars) met[[paste0("Z_", v)]] <- seasonal_z(met[[v]], met$slot)
  met$D_VPD <- met$Z_VPD
  for (v in paste0("SM", 1:4)) met[[paste0("D_", v)]] <- -met[[paste0("Z_", v)]]
  # D_VPD、D_SM均为越大越干；因此正相关表示干湿异常同向变化。
  met$days <- c(pmin(8, as.numeric(diff(met$date))), 8)

  for (v in paste0("SM", 1:4)) {
    met[[paste0("r_", v)]] <- rolling_cor(met$D_VPD, met[[paste0("D_", v)]],
                                         met$date, window_n)
  }
  # 先按相关系数区分“同步耦合”和“弱/反向耦合”；这不代表干旱同时发生。
  met$coupling <- ifelse(is.na(met$r_SM1), NA_character_,
                         ifelse(met$r_SM1 >= cor_cut, "Synchronous", "Asynchronous"))
  # 再单独判断实际干旱是否同时发生，避免将湿润期误称为同步干旱。
  met$drought_state <- NA_character_
  ok <- is.finite(met$D_VPD) & is.finite(met$D_SM1)
  a <- met$D_VPD[ok] >= threshold
  b <- met$D_SM1[ok] >= threshold
  met$drought_state[ok] <- ifelse(a & b, "Both", ifelse(a, "VPD_only",
                                                     ifelse(b, "SM_only", "Neither")))
  met
}

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

drought_analysis_notes <- function(met, gpp, matched, threshold, window_n, cor_cut, bootstrap_n, model_result) {
  two_groups <- model_result$two_groups
  block_n <- model_result$block_n
  notes <- c(
    paste0("实际站点：", unique(met$site), "；题目中的AR-SLu与文件不一致，以数据为准。"),
    sprintf("气象记录%d条，日期%s至%s。", nrow(met), min(met$date), max(met$date)),
    sprintf("GPP记录%d条，精确匹配的有效GPP记录%d条。", nrow(gpp), matched),
    paste0("GPP数据日期：", min(gpp$date), "至", max(gpp$date), "。"),
    "按site和date精确匹配；未插值；source不同是数据来源不同，不作为匹配键。",
    "干旱指数：以全部气象数据为基准，按年内8天时段去季节标准化；SM指数取负号。",
    "本方法是季节性z分数，不是分布拟合SPI/SPEI；GPP异常值的原始基准未知。",
    sprintf("阈值D>=%.1f；持续时间按实际覆盖天数计算，末条最多按8天估计。", threshold),
    "事件end为估计覆盖终日；boundary_or_gap标记可能截断的事件。",
    "强度=时间加权平均D；烈度=sum(D*天数)；同时输出超过阈值部分的累计烈度。",
    sprintf("主分析SM1；向后%d期窗口；r>=%.2f为同步耦合，其余为弱或反向耦合。", window_n, cor_cut),
    if (two_groups) "GPP研究期两类耦合时段均存在，可估计组间差异。" else
      "GPP研究期只有一类耦合时段：无法估计同步/不同步差异，相关结果留空；回归为现有组模型。",
    "相关分组是操作性定义，不是统计显著性判断；干旱同时发生另见drought_state。",
    "回顾性分析使用全期季节基准；不是独立预测验证。",
    "05表按耦合时段比较GPP；06表按实际干旱发生状态比较GPP。",
    "两类均存在时，08表交互项表示同步相对不同步的斜率差；10表缺失组斜率留空。",
    sprintf("回归区间来自%d次、块长%d期移动块bootstrap；不采用普通lm的p值作推断。", bootstrap_n, block_n),
    "区间反映给定指数和分组下的抽样不确定性，不包含基准/分类阈值选择的不确定性。",
    "11表为窗口和土层敏感性描述，不是多重检验后确认的显著结果。",
    "分组比较及回归均为统计关联，不能直接认定干旱导致GPP变化。",
    sprintf("同步组减不同步组GPP均值差=%.4f，95%%块bootstrap区间[%.4f, %.4f]。",
            model_result$mean_difference$GPP_mean_difference,
            model_result$mean_difference$lower_95, model_result$mean_difference$upper_95)
  )
  notes
}

# 完整流程；只导出两份核心CSV，详细统计保留在返回的命名列表中。
# data_dir为输入目录，out_dir为输出目录；相同种子复现bootstrap结果。
# plotter为可选绘图函数；NULL时只计算与导出，标准入口负责传入绘图函数。
run_drought_analysis <- function(data_dir, out_dir, threshold = 1, window_n = 23L,
                                 cor_cut = 0.3, bootstrap_n = 999L,
                                 block_n = 23L, seed = 20261023, plotter = NULL) {
  if (!is.null(plotter) && !is.function(plotter)) stop("plotter必须是函数或NULL")
  set.seed(seed)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  inputs <- read_drought_inputs(data_dir)
  met <- prepare_drought_indices(inputs$met, inputs$vars, threshold, window_n, cor_cut)
  gpp <- inputs$gpp
  indices <- c("D_VPD", paste0("D_SM", 1:4))
  events <- do.call(rbind, lapply(indices, function(v) {
    event_table(met$date, met[[v]], met$days, v, threshold)
  }))

  # 按站点和日期精确连接；GPP缺测保持NA，不插值。
  dat <- merge(met, gpp[c("site", "date", "GPP")],
               by = c("site", "date"), all.x = TRUE)
  dat <- dat[order(dat$date), ]
  matched <- sum(is.finite(dat$GPP))
  if (!matched) stop("两份数据没有日期匹配的有效GPP")
  analysis_data <- dat  # 04表保留原脚本在添加回归协变量前的列结构。
  coupling_summary <- group_summary(dat, "coupling")
  state_summary <- group_summary(dat, "drought_state")
  response <- fit_gpp_response(dat, bootstrap_n, block_n)
  sensitivity <- window_layer_sensitivity(response$data, cor_cut)
  notes <- drought_analysis_notes(met, gpp, matched, threshold, window_n,
                                  cor_cut, bootstrap_n, response)

  stopifnot(nrow(dat) == nrow(met), matched <= nrow(gpp),
            all(abs(met$r_SM1[is.finite(met$r_SM1)]) <= 1 + 1e-12),
            all(events$duration_days > 0),
            all(events$severity_z_days >= events$duration_days * threshold - 1e-10),
            all(response$coefficients$bootstrap_valid >= bootstrap_n * 0.8))

  # 保留分析时间序列和响应斜率，满足研究问题及独立绘图所需。
  save_table(analysis_data, "04_analysis_data.csv", out_dir)
  save_table(response$slopes, "10_GPP_response_slopes.csv", out_dir)
  capture.output(summary(response$model),
                 file = file.path(out_dir, "07_model_descriptive.txt"))
  saveRDS(list(threshold = threshold, cor_cut = cor_cut),
          file.path(out_dir, "plot_settings.rds"))
  if (!is.null(plotter)) {
    plotter(dat = response$data, slopes = response$slopes,
            threshold = threshold, cor_cut = cor_cut,
            two_groups = response$two_groups, out_dir = file.path(out_dir, "images"))
  }
  writeLines(notes, file.path(out_dir, "15_methods_and_results.txt"), useBytes = TRUE)
  capture.output(sessionInfo(), file = file.path(out_dir, "16_sessionInfo.txt"))

  cat("运行完成！结果目录：", out_dir, "\n")
  print(coupling_summary)
  print(state_summary)
  print(response$slopes)
  invisible(list(data = response$data, events = events,
                 coupling_summary = coupling_summary, state_summary = state_summary,
                 model = response$model, coefficients = response$coefficients,
                 mean_difference = response$mean_difference, slopes = response$slopes,
                 sensitivity = sensitivity, quality = inputs$quality, notes = notes))
}
