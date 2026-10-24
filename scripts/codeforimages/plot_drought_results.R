# 绘图函数库：加载时只定义函数，不读取数据、不输出图片。

# 从完整分析的输出恢复绘图数据，不重新执行分析。
read_plot_inputs <- function(out_dir) {
  files <- c("04_analysis_data.csv", file.path("metadata", "plot_settings.rds"))
  missing <- files[!file.exists(file.path(out_dir, files))]
  if (length(missing)) {
    stop("缺少绘图输入：", paste(missing, collapse = ", "),
         "。请先运行scripts/run_drought_analysis.R。", call. = FALSE)
  }
  dat <- read.csv(file.path(out_dir, files[1L]), stringsAsFactors = FALSE,
                  na.strings = c("", "NA"))
  settings <- readRDS(file.path(out_dir, files[2L]))
  required_data <- c("date", "site", "D_VPD", "D_SM1", "r_SM1", "GPP",
                     "coupling", "drought_state")
  if (!all(required_data %in% names(dat))) stop("绘图CSV缺少必要列，请重新运行完整分析。")
  dat$date <- as.Date(dat$date)
  if (anyNA(dat$date)) stop("绘图数据包含无效日期。", call. = FALSE)
  dat <- dat[order(dat$date), ]
  dat$coupling <- factor(dat$coupling, levels = c("Asynchronous", "Synchronous"))
  groups <- unique(as.character(dat$coupling[is.finite(dat$GPP) & !is.na(dat$coupling)]))
  list(dat = dat, threshold = settings$threshold,
       cor_cut = settings$cor_cut, two_groups = length(groups) == 2L)
}

# 英文图标签保证跨系统字体可用，分析注释与方法说明为中文。
draw_plot <- function(name, fun, out_dir, height = 1000) {
  png(file.path(out_dir, name), width = 1600, height = height, res = 150)
  on.exit(dev.off())
  fun()
}

plot_drought_results <- function(dat, threshold, cor_cut, two_groups, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dat$coupling <- factor(dat$coupling, levels = c("Asynchronous", "Synchronous"))
  draw_plot(out_dir = out_dir, name = "12_time_series.png", function() {
    par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
    plot(dat$date, dat$D_VPD, type = "l", col = "firebrick", ylab = "Drought z-score",
         xlab = "", main = paste(unique(dat$site), "- VPD and surface SM"),
         ylim = range(c(dat$D_VPD, dat$D_SM1), na.rm = TRUE))
    lines(dat$date, dat$D_SM1, col = "steelblue")
    abline(h = threshold, lty = 2)
    legend("topleft", c("VPD (higher = drier)", "SM1 (higher = drier)"),
           col = c("firebrick", "steelblue"), lty = 1, bty = "n", cex = 0.8)
    plot(dat$date, dat$r_SM1, type = "l", ylim = c(-1, 1), ylab = "Rolling r", xlab = "")
    abline(h = c(0, cor_cut), lty = c(3, 2), col = c("grey50", "firebrick"))
    plot(dat$date, dat$GPP, type = "l", ylab = "GPP anomaly", xlab = "Date")
    abline(h = 0, lty = 2)
  }, height = 1200)
  draw_plot(out_dir = out_dir, name = "13_GPP_groups.png", function() {
    par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
    boxplot(GPP ~ coupling, data = dat, col = c("grey80", "tan"), las = 2,
            xlab = "", ylab = "GPP anomaly", main = "By rolling-correlation regime")
    if (!two_groups) mtext("Only one regime observed; no between-regime comparison", side = 3, cex = 0.65)
    abline(h = 0, lty = 2)
    boxplot(GPP ~ drought_state, data = dat,
            col = c("tomato", "grey85", "steelblue", "gold"), las = 2,
            xlab = "", ylab = "GPP anomaly", main = "By actual drought occurrence")
    abline(h = 0, lty = 2)
  })
}
