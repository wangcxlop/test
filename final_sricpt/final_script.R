# 从项目根目录运行：Rscript final_sricpt/final_script.R
# 首次使用：install.packages("dplyr")
suppressPackageStartupMessages(library(dplyr))
data_dir <- "Rbackage"
out_dir <- "output"
threshold <- 1
window_n <- 23L
cor_cut <- 0.3

# 标准化；至少5个有效样本且标准差非零
z_score <- function(x) {
  z <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) >= 5 && sd(x[ok]) > 0) z[ok] <- (x[ok] - mean(x[ok])) / sd(x[ok])
  z
}

# 相关系数；至少3个有效配对且变量不恒定
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}

# 1. 读取数据
met <- read.csv(file.path(data_dir, "ERA5L_TWSA.csv")) |>
  select(site, date, VPD, SM1) |> mutate(date = as.Date(date)) |> arrange(site, date)
gpp <- read.csv(file.path(data_dir, "GPP_anomaly.csv")) |>
  select(site, date, GPP) |> mutate(date = as.Date(date)) |> arrange(site, date)
if (anyNA(met$date) || anyNA(gpp$date) || anyNA(met$site) || anyNA(gpp$site)) stop("日期或站点无效")
if (n_distinct(met$site) != 1 || n_distinct(gpp$site) != 1 || unique(met$site) != unique(gpp$site)) stop("需要同一站点的数据")
if (any(count(met, site, date)$n > 1) || any(count(gpp, site, date)$n > 1)) stop("站点日期有重复")
if (!is.numeric(met$VPD) || !is.numeric(met$SM1) || !is.numeric(gpp$GPP)) stop("分析变量必须为数值")

# 2. 按年内8天时段标准化，两个指数均为越大越干
met <- met |>
  mutate(slot = (as.integer(format(date, "%j")) - 1L) %/% 8L + 1L) |>
  group_by(slot) |> mutate(Z_VPD = z_score(VPD), Z_SM1 = z_score(SM1)) |> ungroup() |>
  arrange(date) |>
  mutate(D_VPD = Z_VPD, D_SM1 = -Z_SM1, days = c(pmin(8, as.numeric(diff(date))), 8)) |>
  as.data.frame()

# 3. 统计连续干旱的强度、时长和烈度
indices <- c("D_VPD", "D_SM1")
events <- data.frame(index = character(), start = as.Date(character()), end = as.Date(character()),
                     n_periods = integer(), duration_days = numeric(), intensity = numeric(),
                     severity_z_days = numeric(), boundary_or_gap = logical())
gap <- c(TRUE, as.numeric(diff(met$date)) > 8)
for (index in indices) {
  x <- met[[index]]
  dry <- is.finite(x) & x >= threshold
  starts <- which(dry & (gap | !c(FALSE, head(dry, -1))))
  for (i in starts) {
    j <- i
    while (j < nrow(met) && dry[j + 1] && !gap[j + 1]) j <- j + 1
    ids <- i:j
    boundary <- i == 1 || j == nrow(met)
    if (i > 1) boundary <- boundary || gap[i] || !is.finite(x[i - 1])
    if (j < nrow(met)) boundary <- boundary || gap[j + 1] || !is.finite(x[j + 1])
    events <- bind_rows(events, data.frame(
      index = index, start = met$date[i], end = met$date[j] + met$days[j] - 1,
      n_periods = length(ids), duration_days = sum(met$days[ids]),
      intensity = weighted.mean(x[ids], met$days[ids]),
      severity_z_days = sum(x[ids] * met$days[ids]), boundary_or_gap = boundary))
  }
}
event_summary <- events |> mutate(index = factor(index, levels = indices)) |>
  group_by(index, .drop = FALSE) |>
  summarise(n_events = n(),
            mean_intensity = if (n() > 0) mean(intensity) else NA_real_,
            mean_duration_days = if (n() > 0) mean(duration_days) else NA_real_,
            mean_severity_z_days = if (n() > 0) mean(severity_z_days) else NA_real_, .groups = "drop") |>
  mutate(index = as.character(index)) |> as.data.frame()

# 4. 向后滑动相关；至少80%有效配对，不跨越8天以上缺口
met$r_SM1 <- NA_real_
if (nrow(met) >= window_n) {
  for (i in seq.int(window_n, nrow(met))) {
    ids <- seq.int(i - window_n + 1, i)
    if (any(as.numeric(diff(met$date[ids])) > 8)) next
    x <- met$D_VPD[ids]
    y <- met$D_SM1[ids]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) >= ceiling(window_n * 0.8) && sd(x[ok]) > 0 && sd(y[ok]) > 0) {
      met$r_SM1[i] <- cor(x[ok], y[ok])
    }
  }
}
# 耦合同步不等于同时干旱
met <- met |> mutate(
  coupling = case_when(is.na(r_SM1) ~ NA_character_, r_SM1 >= cor_cut ~ "Synchronous", TRUE ~ "Asynchronous"),
  drought_state = case_when(
    !is.finite(D_VPD) | !is.finite(D_SM1) ~ NA_character_,
    D_VPD >= threshold & D_SM1 >= threshold ~ "Both",
    D_VPD >= threshold ~ "VPD_only", D_SM1 >= threshold ~ "SM_only", TRUE ~ "Neither"))

# 5. 合并GPP并计算分组均值与相关
# GPP已经是异常值，不再标准化
dat <- met |> left_join(gpp, by = c("site", "date")) |> arrange(date) |>
  mutate(threshold = .env$threshold, cor_cut = .env$cor_cut)
if (!any(is.finite(dat$GPP))) stop("没有日期匹配的有效GPP")
levels_by_group <- list(coupling = c("Asynchronous", "Synchronous"),
                        drought_state = c("Both", "VPD_only", "SM_only", "Neither"))
response <- data.frame()
for (kind in names(levels_by_group)) {
  result <- dat |> rename(group = all_of(kind)) |>
    mutate(group = factor(group, levels = levels_by_group[[kind]])) |>
    filter(!is.na(group), is.finite(GPP)) |> group_by(group, .drop = FALSE) |>
    summarise(n_GPP = n(), mean_GPP = if (n() > 0) mean(GPP) else NA_real_,
              n_VPD = sum(is.finite(D_VPD)), r_VPD = safe_cor(GPP, D_VPD),
              n_SM1 = sum(is.finite(D_SM1)), r_SM1 = safe_cor(GPP, D_SM1), .groups = "drop") |>
    mutate(grouping = kind, group = as.character(group)) |> relocate(grouping, group) |> as.data.frame()
  response <- bind_rows(response, result)
}

# 6. 保存表格和简短说明
dir.create(file.path(out_dir, "reports"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "images"), recursive = TRUE, showWarnings = FALSE)
write.csv(dat, file.path(out_dir, "04_analysis_data.csv"), row.names = FALSE, na = "")
write.csv(response, file.path(out_dir, "10_GPP_group_summary.csv"), row.names = FALSE, na = "")
notes <- capture.output({
  cat("站点：", unique(dat$site), "；按站点和日期精确匹配，不插值。\n")
  cat("VPD按年内8天时段标准化，SM1标准化后取负号；均为越大越干。\n")
  cat("干旱阈值：", threshold, "；滑动窗口：", window_n, "期；同步耦合阈值：", cor_cut, "。\n")
  cat("强度为时间加权平均指数，时长为覆盖天数，烈度为sum(指数*天数)。\n")
  cat("\n干旱事件汇总（各事件等权平均）：\n")
  print(event_summary, row.names = FALSE, digits = 4)
  cat("\n分组GPP均值与Pearson相关：\n")
  print(response, row.names = FALSE, digits = 4)
  missing <- response$group[response$grouping == "coupling" & response$n_GPP == 0]
  if (length(missing)) cat("\n没有有效GPP样本的耦合组：", missing, "；无法比较同步与不同步组。\n")
  cat("NA表示无法估计。耦合同步不等于同时发生干旱；这些结果是关联，不代表因果或统计显著。\n")
})
writeLines(notes, file.path(out_dir, "reports", "15_methods_and_results.txt"), useBytes = TRUE)

# 7. 绘制两张图
plot_data <- dat |> mutate(coupling = factor(coupling, levels = c("Asynchronous", "Synchronous")))
two_groups <- sum(response$grouping == "coupling" & response$n_GPP > 0) == 2
png(file.path(out_dir, "images", "12_time_series.png"), width = 1600, height = 1200, res = 150)
par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
plot(plot_data$date, plot_data$D_VPD, type = "l", col = "firebrick", ylab = "Drought z-score",
     xlab = "", main = paste(unique(plot_data$site), "- VPD and surface SM"),
     ylim = range(c(plot_data$D_VPD, plot_data$D_SM1), na.rm = TRUE))
lines(plot_data$date, plot_data$D_SM1, col = "steelblue")
abline(h = threshold, lty = 2)
legend("topleft", c("VPD (higher = drier)", "SM1 (higher = drier)"),
       col = c("firebrick", "steelblue"), lty = 1, bty = "n", cex = 0.8)
plot(plot_data$date, plot_data$r_SM1, type = "l", ylim = c(-1, 1), ylab = "Rolling r", xlab = "")
abline(h = c(0, cor_cut), lty = c(3, 2), col = c("grey50", "firebrick"))
plot(plot_data$date, plot_data$GPP, type = "l", ylab = "GPP anomaly", xlab = "Date")
abline(h = 0, lty = 2)
dev.off()
png(file.path(out_dir, "images", "13_GPP_groups.png"), width = 1600, height = 1000, res = 150)
par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
boxplot(GPP ~ coupling, data = plot_data, col = c("grey80", "tan"), las = 2,
        xlab = "", ylab = "GPP anomaly", main = "By rolling-correlation regime")
if (!two_groups) mtext("Only one regime observed; no between-regime comparison", side = 3, cex = 0.65)
abline(h = 0, lty = 2)
boxplot(GPP ~ drought_state, data = plot_data, col = c("tomato", "grey85", "steelblue", "gold"),
        las = 2, xlab = "", ylab = "GPP anomaly", main = "By actual drought occurrence")
abline(h = 0, lty = 2)
dev.off()
print(event_summary)
print(response)
cat("完成！结果保存在", out_dir, "\n")
