# 读取数据、计算指标、分组比较并保存
run_drought_analysis <- function(data_dir, out_dir, threshold = 1, window_n = 23L, cor_cut = 0.3) {
  inputs <- read_drought_inputs(data_dir)
  met <- prepare_drought_indices(inputs$met, threshold, window_n, cor_cut)
  events <- data.frame()
  event_summary <- data.frame()
  for (index in c("D_VPD", "D_SM1")) {
    e <- event_table(met$date, met[[index]], met$days, index, threshold)
    events <- rbind(events, e)
    event_summary <- rbind(event_summary, data.frame(
      index = index, n_events = nrow(e),
      mean_intensity = if (nrow(e)) mean(e$intensity) else NA_real_,
      mean_duration_days = if (nrow(e)) mean(e$duration_days) else NA_real_,
      mean_severity_z_days = if (nrow(e)) mean(e$severity_z_days) else NA_real_))
  }

  dat <- merge(met, inputs$gpp, by = c("site", "date"), all.x = TRUE)
  dat <- dat[order(dat$date), ]
  if (!any(is.finite(dat$GPP))) stop("两份数据没有日期匹配的有效GPP")
  dat$threshold <- threshold
  dat$cor_cut <- cor_cut
  response <- analyze_gpp_response(dat)
  write_drought_outputs(dat, response, event_summary, out_dir, window_n)
  print(event_summary)
  print(response)
  cat("分析完成，结果保存在：", out_dir, "\n")
  invisible(list(data = dat, events = events, group_summary = response, event_summary = event_summary))
}
