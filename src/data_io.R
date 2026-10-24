# 数据读写与质量检查；加载时只定义函数。
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
  vars <- c("VPD", "SM1")
  met <- read_input("ERA5L_TWSA.csv", c("site", "date", vars), data_dir)
  gpp <- read_input("GPP_anomaly.csv", c("site", "date", "GPP"), data_dir)
  if (length(unique(met$site)) != 1L || length(unique(gpp$site)) != 1L ||
      unique(met$site) != unique(gpp$site)) stop("本脚本要求两份文件为同一个站点")
  if (!all(vapply(met[vars], is.numeric, logical(1))) || !is.numeric(gpp$GPP)) {
    stop("分析变量必须为数值列")
  }
  list(met = met[c("site", "date", vars)], gpp = gpp[c("site", "date", "GPP")], vars = vars)
}

# 导出两份核心CSV，报告和元数据分别写入子目录。
# 统计计算在其他模块完成，这里只负责文件写入。
write_drought_outputs <- function(analysis_data, response, notes, threshold, cor_cut, out_dir) {
  report_dir <- file.path(out_dir, "reports")
  metadata_dir <- file.path(out_dir, "metadata")
  for (path in c(out_dir, report_dir, metadata_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  save_table(analysis_data, "04_analysis_data.csv", out_dir)
  save_table(response, "10_GPP_group_summary.csv", out_dir)
  writeLines(notes, file.path(report_dir, "15_methods_and_results.txt"), useBytes = TRUE)
  saveRDS(list(threshold = threshold, cor_cut = cor_cut),
          file.path(metadata_dir, "plot_settings.rds"))
}
