# 绘制时间序列和GPP箱线图
plot_drought_results <- function(dat, out_dir) {
  threshold <- dat$threshold[1]
  cor_cut <- dat$cor_cut[1]
  groups <- unique(dat$coupling[is.finite(dat$GPP) & !is.na(dat$coupling)])
  two_groups <- length(groups) == 2L
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dat$coupling <- factor(dat$coupling, levels = c("Asynchronous", "Synchronous"))
  png(file.path(out_dir, "12_time_series.png"), width = 1600, height = 1200, res = 150)
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
  dev.off()
  png(file.path(out_dir, "13_GPP_groups.png"), width = 1600, height = 1000, res = 150)
  par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
  boxplot(GPP ~ coupling, data = dat, col = c("grey80", "tan"), las = 2,
          xlab = "", ylab = "GPP anomaly", main = "By rolling-correlation regime")
  if (!two_groups) mtext("Only one regime observed; no between-regime comparison", side = 3, cex = 0.65)
  abline(h = 0, lty = 2)
  boxplot(GPP ~ drought_state, data = dat,
          col = c("tomato", "grey85", "steelblue", "gold"), las = 2,
          xlab = "", ylab = "GPP anomaly", main = "By actual drought occurrence")
  abline(h = 0, lty = 2)
  dev.off()
}
