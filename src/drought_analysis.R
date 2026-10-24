# 读取数据、计算指标、分组比较并保存
run_drought_analysis <- function(data_dir, out_dir, threshold = 1, window_n = 23L, cor_cut = 0.3) {
  inputs <- read_drought_inputs(data_dir)
  met <- prepare_drought_indices(inputs$met, threshold, window_n, cor_cut)
  events <- list()
  for (index in c("D_VPD", "D_SM1")) {
    events[[index]] <- event_table(met$date, met[[index]], met$days, index, threshold)
  }
  events <- bind_rows(events)
  event_summary <- events |>
    mutate(index = factor(index, levels = c("D_VPD", "D_SM1"))) |>
    group_by(index, .drop = FALSE) |>
    summarise(
      n_events = n(),
      mean_intensity = if (n() > 0) mean(intensity) else NA_real_,
      mean_duration_days = if (n() > 0) mean(duration_days) else NA_real_,
      mean_severity_z_days = if (n() > 0) mean(severity_z_days) else NA_real_,
      .groups = "drop") |>
    mutate(index = as.character(index)) |>
    as.data.frame()

  dat <- met |>
    left_join(inputs$gpp, by = c("site", "date")) |>
    arrange(date) |>
    mutate(threshold = .env$threshold, cor_cut = .env$cor_cut)
  if (!any(is.finite(dat$GPP))) stop("两份数据没有日期匹配的有效GPP")
  response <- analyze_gpp_response(dat)
  write_drought_outputs(dat, response, event_summary, out_dir, window_n)
  print(event_summary)
  print(response)
  cat("分析完成，结果保存在：", out_dir, "\n")
  invisible(list(data = dat, events = events, group_summary = response, event_summary = event_summary))
}
