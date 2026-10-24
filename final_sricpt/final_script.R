# 项目根目录运行：Rscript final_sricpt/final_script.R
suppressPackageStartupMessages(library(dplyr))
threshold <- 1
window_n <- 23L
cor_cut <- 0.3
out_dir <- "output"

# 季节标准化与有效配对相关
z_score <- function(x) {
  z <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) >= 5 && sd(x[ok]) > 0) z[ok] <- (x[ok] - mean(x[ok])) / sd(x[ok])
  z
}
safe_cor <- function(x, y, min_pairs = 3) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_pairs || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}

# 读取并标准化VPD和SM1，数值越大越干
met <- read.csv("Rbackage/ERA5L_TWSA.csv") |>
  select(site, date, VPD, SM1) |> mutate(date = as.Date(date)) |> arrange(date)
gpp <- read.csv("Rbackage/GPP_anomaly.csv") |> select(site, date, GPP) |> mutate(date = as.Date(date))
met <- met |> mutate(slot = (as.integer(format(date, "%j")) - 1L) %/% 8L + 1L) |>
  group_by(slot) |> mutate(D_VPD = z_score(VPD), D_SM1 = -z_score(SM1)) |> ungroup() |>
  arrange(date) |> mutate(days = c(pmin(8, as.numeric(diff(date))), 8))
gap <- c(TRUE, as.numeric(diff(met$date)) > 8)

# 连续干旱的时长、强度和烈度
event_tables <- list()
for (v in c("D_VPD", "D_SM1")) {
  event_tables[[v]] <- met |>
    mutate(index = v, D = .data[[v]], dry = is.finite(D) & D >= threshold, event = cumsum(!dry | gap)) |>
    filter(dry) |> group_by(index, event) |>
    summarise(start = first(date), end = last(date) + last(days) - 1,
              duration_days = sum(days), intensity = weighted.mean(D, days),
              severity_z_days = sum(D * days), .groups = "drop") |> select(-event)
}
events <- bind_rows(event_tables)

# 滑动相关：80%有效配对，不跨越8天以上缺口
met$r_SM1 <- NA_real_
for (i in seq_len(max(0, nrow(met) - window_n + 1))) {
  ids <- seq.int(i, length.out = window_n)
  if (any(as.numeric(diff(met$date[ids])) > 8)) next
  met$r_SM1[i + window_n - 1] <- safe_cor(met$D_VPD[ids], met$D_SM1[ids], max(2, ceiling(window_n * 0.8)))
}
dat <- met |> left_join(gpp, by = c("site", "date")) |>
  mutate(coupling = factor(if_else(r_SM1 >= cor_cut, "Synchronous", "Asynchronous"),
                           levels = c("Asynchronous", "Synchronous"))) |>
  select(site, date, D_VPD, D_SM1, r_SM1, coupling, GPP)
if (!any(is.finite(dat$GPP))) stop("没有日期匹配的有效GPP")
response <- dat |> filter(!is.na(coupling), is.finite(GPP)) |> group_by(coupling, .drop = FALSE) |>
  summarise(n_GPP = n(), mean_GPP = if (n() > 0) mean(GPP) else NA_real_,
            r_VPD = safe_cor(GPP, D_VPD), r_SM1 = safe_cor(GPP, D_SM1), .groups = "drop")

# 保存结果；空组不能比较，相关不代表因果
dir.create(file.path(out_dir, "reports"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "images"), recursive = TRUE, showWarnings = FALSE)
write.csv(dat, file.path(out_dir, "04_analysis_data.csv"), row.names = FALSE, na = "")
write.csv(response, file.path(out_dir, "10_GPP_group_summary.csv"), row.names = FALSE, na = "")
capture.output({
  cat("干旱阈值：", threshold, "；窗口：", window_n, "；同步阈值：", cor_cut, "\n", sep = "")
  print(as.data.frame(events), row.names = FALSE)
  print(as.data.frame(response), row.names = FALSE)
  if (any(response$n_GPP == 0)) cat("存在无有效GPP样本的组，无法比较同步与不同步组。\n")
  cat("耦合同步不等于同时干旱。NA表示无法估计，相关不代表因果或统计显著。\n")
}, file = file.path(out_dir, "reports", "15_methods_and_results.txt"))

# 三幅时间序列和一幅GPP分组箱线图
png(file.path(out_dir, "images", "12_time_series.png"), width = 1600, height = 1200, res = 150)
par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
plot(dat$date, dat$D_VPD, type = "l", col = "red", xlab = "", ylab = "Drought index",
     ylim = range(c(-1, 1, dat$D_VPD, dat$D_SM1), na.rm = TRUE))
lines(dat$date, dat$D_SM1, col = "blue")
legend("topleft", c("VPD", "SM1"), col = c("red", "blue"), lty = 1, bty = "n")
abline(h = threshold, lty = 2)
plot(dat$date, dat$r_SM1, type = "l", ylim = c(-1, 1), xlab = "", ylab = "Rolling r")
abline(h = cor_cut, lty = 2)
plot(dat$date, dat$GPP, type = "l", xlab = "Date", ylab = "GPP anomaly")
dev.off()
png(file.path(out_dir, "images", "13_GPP_groups.png"), width = 1400, height = 900, res = 150)
par(mar = c(7, 4, 3, 1))
if (any(response$n_GPP > 0)) boxplot(GPP ~ coupling, data = filter(dat, is.finite(GPP)), las = 2, xlab = "", ylab = "GPP anomaly") else plot.new()
if (any(response$n_GPP == 0)) mtext("No valid GPP in one or both groups", side = 3, cex = 0.8)
dev.off()
print(response)
