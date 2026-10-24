# 读取和保存数据
# GPP已是异常值，无需再去季节化

# 检查日期、站点和重复记录
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
  list(met = met[c("site", "date", vars)], gpp = gpp[c("site", "date", "GPP")])
}

# 保存结果
write_drought_outputs <- function(dat, response, event_summary, out_dir, window_n) {
  dir.create(file.path(out_dir, "reports"), recursive = TRUE, showWarnings = FALSE)
  write.csv(dat, file.path(out_dir, "04_analysis_data.csv"), row.names = FALSE, na = "")
  write.csv(response, file.path(out_dir, "10_GPP_group_summary.csv"), row.names = FALSE, na = "")
  notes <- capture.output({
    cat("站点：", unique(dat$site), "；按站点和日期精确匹配，不插值。\n")
    cat("VPD按年内8天时段标准化，SM1标准化后取负号；均为越大越干。\n")
    cat("干旱阈值：", dat$threshold[1], "；滑动窗口：", window_n,
        "期；同步耦合阈值：", dat$cor_cut[1], "。\n")
    cat("强度为时间加权平均指数，时长为覆盖天数，烈度为sum(指数*天数)。\n")
    cat("\n干旱事件汇总（各事件等权平均）：\n")
    print(event_summary, row.names = FALSE, digits = 4)
    cat("\n分组GPP均值与Pearson相关：\n")
    print(response, row.names = FALSE, digits = 4)
    missing <- response$group[response$grouping == "coupling" & response$n_GPP == 0]
    if (length(missing)) cat("\n没有有效GPP样本的耦合组：", missing,
                             "；无法比较同步与不同步组。\n")
    cat("NA表示无法估计。耦合同步不等于同时发生干旱；这些结果是关联，不代表因果或统计显著。\n")
  })
  writeLines(notes, file.path(out_dir, "reports", "15_methods_and_results.txt"), useBytes = TRUE)
}
