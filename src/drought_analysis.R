# 干旱分析主流程与中文说明；加载时只定义函数。
# 完整分析入口负责加载四个分析模块和绘图模块。

drought_analysis_notes <- function(met, gpp, matched, threshold, window_n, cor_cut, events, response) {
  event_summary <- do.call(rbind, lapply(c("D_VPD", "D_SM1"), function(v) {
    e <- events[events$index == v, ]
    data.frame(index = v, n_events = nrow(e),
               mean_intensity = if (nrow(e)) mean(e$intensity) else NA_real_,
               mean_duration_days = if (nrow(e)) mean(e$duration_days) else NA_real_,
               mean_severity_z_days = if (nrow(e)) mean(e$severity_z_days) else NA_real_)
  }))
  two_groups <- sum(response$grouping == "coupling" & response$n_GPP > 0) == 2L
  missing_groups <- response$group[response$grouping == "coupling" & response$n_GPP == 0]
  missing_labels <- c(Asynchronous = "不同步组", Synchronous = "同步组")[missing_groups]
  both <- response[response$grouping == "drought_state" & response$group == "Both", ]
  neither <- response[response$grouping == "drought_state" & response$group == "Neither", ]
  occurrence_result <- if (both$n_GPP > 0 && neither$n_GPP > 0) {
    sprintf("共同干旱与均无干旱时段的GPP均值分别为%.4f和%.4f，均值差为%.4f。",
            both$mean_GPP, neither$mean_GPP, both$mean_GPP - neither$mean_GPP)
  } else "共同干旱或均无干旱组缺少有效GPP，不能比较这两组均值。"
  c(
    paste0("站点：", unique(met$site), "；题目站点名称与文件不一致，以实际数据为准。"),
    sprintf("气象数据：%s至%s，共%d条。", min(met$date), max(met$date), nrow(met)),
    sprintf("GPP数据：%s至%s，共%d条；有效日期匹配%d条。",
            min(gpp$date), max(gpp$date), nrow(gpp), matched),
    "按站点和日期精确匹配，不插值；GPP使用原始异常值，其计算基准未知。",
    "VPD和SM1按年内8天时段、以全期多年均值和标准差计算季节性z分数。",
    "D_VPD=z(VPD)，D_SM1=-z(SM1)，均为越大越干；这不是分布拟合SPI/SPEI。",
    sprintf("D>=%.1f识别干旱；强度为事件内指数的时间加权均值，持续时间为实际覆盖天数，烈度为sum(D*天数)。", threshold),
    "末条记录覆盖时长按8天估计；boundary_or_gap标记可能截断的事件；事件详情见返回对象events。",
    sprintf("向后%d期滑动相关，至少80%%有效配对，不跨越超过8天的日期缺口。", window_n),
    sprintf("r>=%.2f定义同步耦合，其余为弱/反向耦合；这是操作性分类，不是显著性检验。", cor_cut),
    "coupling描述干湿异常的耦合；drought_state另外区分Both（共同干旱）、VPD_only、SM_only和Neither（均无干旱）。",
    "分组统计：有效GPP样本数、GPP均值、GPP与两类干旱指数的Pearson相关及有效配对数。",
    "相关系数至少需要3个有效配对且两个变量均有变化；NA表示无法估计，不代表零相关。",
    "", "干旱事件汇总（各事件等权平均）：",
    capture.output(print(event_summary, row.names = FALSE, digits = 4)),
    "", "分组GPP与相关结果：",
    capture.output(print(response, row.names = FALSE, digits = 4)),
    "", "结论与限制：",
    if (two_groups) "两类耦合时段均有有效GPP，均值及相关差异见上表。" else
      paste0("当前规则下", paste(missing_labels, collapse = "、"),
             "没有有效GPP样本，无法比较同步与不同步组；不改变阈值制造分组。"),
    occurrence_result,
    "相关系数表示组内干旱指数与GPP的同期关联；相关和均值差异不等于因果影响，也不表示统计显著。"
  )
}

run_drought_analysis <- function(data_dir, out_dir, threshold = 1, window_n = 23L,
                                 cor_cut = 0.3, plotter = NULL) {
  if (!is.null(plotter) && !is.function(plotter)) stop("plotter必须是函数或NULL")
  inputs <- read_drought_inputs(data_dir)
  met <- prepare_drought_indices(inputs$met, inputs$vars, threshold, window_n, cor_cut)
  gpp <- inputs$gpp
  events <- do.call(rbind, lapply(c("D_VPD", "D_SM1"), function(v) {
    event_table(met$date, met[[v]], met$days, v, threshold)
  }))
  dat <- merge(met, gpp, by = c("site", "date"), all.x = TRUE)
  dat <- dat[order(dat$date), ]
  matched <- sum(is.finite(dat$GPP))
  if (!matched) stop("两份数据没有日期匹配的有效GPP")
  response <- analyze_gpp_response(dat)
  notes <- drought_analysis_notes(met, gpp, matched, threshold, window_n, cor_cut, events, response)
  write_drought_outputs(dat, response, notes, threshold, cor_cut, out_dir)
  if (!is.null(plotter)) {
    plotter(dat = dat, threshold = threshold, cor_cut = cor_cut,
            two_groups = sum(response$grouping == "coupling" & response$n_GPP > 0) == 2L,
            out_dir = file.path(out_dir, "images"))
  }
  cat("运行完成！结果目录：", out_dir, "\n")
  print(response)
  invisible(list(data = dat, events = events, group_summary = response, notes = notes))
}
