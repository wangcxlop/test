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

# 导出两份核心CSV，报告和元数据分别写入子目录。
# 统计计算在其他模块完成，这里只负责文件写入。
write_drought_outputs <- function(analysis_data, response, notes, threshold, cor_cut, out_dir) {
  report_dir <- file.path(out_dir, "reports")
  metadata_dir <- file.path(out_dir, "metadata")
  for (path in c(out_dir, report_dir, metadata_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  save_table(analysis_data, "04_analysis_data.csv", out_dir)
  save_table(response$slopes, "10_GPP_response_slopes.csv", out_dir)
  capture.output(summary(response$model),
                 file = file.path(report_dir, "07_model_descriptive.txt"))
  writeLines(notes, file.path(report_dir, "15_methods_and_results.txt"), useBytes = TRUE)
  saveRDS(list(threshold = threshold, cor_cut = cor_cut),
          file.path(metadata_dir, "plot_settings.rds"))
  capture.output(sessionInfo(), file = file.path(metadata_dir, "16_sessionInfo.txt"))
}
