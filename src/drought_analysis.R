# 干旱分析主流程与中文方法说明；加载时只定义函数。
# 运行前依次加载data_io.R、drought_indices.R、gpp_response.R和本文件。
# 日常运行请使用scripts/run_drought_analysis.R，它会按顺序加载全部模块。

drought_analysis_notes <- function(met, gpp, matched, threshold, window_n, cor_cut, bootstrap_n, model_result) {
  two_groups <- model_result$two_groups
  block_n <- model_result$block_n
  notes <- c(
    paste0("实际站点：", unique(met$site), "；题目中的AR-SLu与文件不一致，以数据为准。"),
    sprintf("气象记录%d条，日期%s至%s。", nrow(met), min(met$date), max(met$date)),
    sprintf("GPP记录%d条，精确匹配的有效GPP记录%d条。", nrow(gpp), matched),
    paste0("GPP数据日期：", min(gpp$date), "至", max(gpp$date), "。"),
    "按site和date精确匹配；未插值；source不同是数据来源不同，不作为匹配键。",
    "干旱指数：以全部气象数据为基准，按年内8天时段去季节标准化；SM指数取负号。",
    "本方法是季节性z分数，不是分布拟合SPI/SPEI；GPP异常值的原始基准未知。",
    sprintf("阈值D>=%.1f；持续时间按实际覆盖天数计算，末条最多按8天估计。", threshold),
    "事件end为估计覆盖终日；boundary_or_gap标记可能截断的事件。",
    "强度=时间加权平均D；烈度=sum(D*天数)；同时输出超过阈值部分的累计烈度。",
    sprintf("主分析SM1；向后%d期窗口；r>=%.2f为同步耦合，其余为弱或反向耦合。", window_n, cor_cut),
    if (two_groups) "GPP研究期两类耦合时段均存在，可估计组间差异。" else
      "GPP研究期只有一类耦合时段：无法估计同步/不同步差异，相关结果留空；回归为现有组模型。",
    "相关分组是操作性定义，不是统计显著性判断；干旱同时发生另见drought_state。",
    "回顾性分析使用全期季节基准；不是独立预测验证。",
    "05表按耦合时段比较GPP；06表按实际干旱发生状态比较GPP。",
    "两类均存在时，08表交互项表示同步相对不同步的斜率差；10表缺失组斜率留空。",
    sprintf("回归区间来自%d次、块长%d期移动块bootstrap；不采用普通lm的p值作推断。", bootstrap_n, block_n),
    "区间反映给定指数和分组下的抽样不确定性，不包含基准/分类阈值选择的不确定性。",
    "11表为窗口和土层敏感性描述，不是多重检验后确认的显著结果。",
    "分组比较及回归均为统计关联，不能直接认定干旱导致GPP变化。",
    sprintf("同步组减不同步组GPP均值差=%.4f，95%%块bootstrap区间[%.4f, %.4f]。",
            model_result$mean_difference$GPP_mean_difference,
            model_result$mean_difference$lower_95, model_result$mean_difference$upper_95)
  )
  notes
}

# 完整流程；只导出两份核心CSV，详细统计保留在返回的命名列表中。
# data_dir为输入目录，out_dir为输出目录；相同种子复现bootstrap结果。
# plotter为可选绘图函数；NULL时只计算与导出，标准入口负责传入绘图函数。
run_drought_analysis <- function(data_dir, out_dir, threshold = 1, window_n = 23L,
                                 cor_cut = 0.3, bootstrap_n = 999L,
                                 block_n = 23L, seed = 20261023, plotter = NULL) {
  if (!is.null(plotter) && !is.function(plotter)) stop("plotter必须是函数或NULL")
  set.seed(seed)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  inputs <- read_drought_inputs(data_dir)
  met <- prepare_drought_indices(inputs$met, inputs$vars, threshold, window_n, cor_cut)
  gpp <- inputs$gpp
  indices <- c("D_VPD", paste0("D_SM", 1:4))
  events <- do.call(rbind, lapply(indices, function(v) {
    event_table(met$date, met[[v]], met$days, v, threshold)
  }))

  # 按站点和日期精确连接；GPP缺测保持NA，不插值。
  dat <- merge(met, gpp[c("site", "date", "GPP")],
               by = c("site", "date"), all.x = TRUE)
  dat <- dat[order(dat$date), ]
  matched <- sum(is.finite(dat$GPP))
  if (!matched) stop("两份数据没有日期匹配的有效GPP")
  analysis_data <- dat  # 04表保留原脚本在添加回归协变量前的列结构。
  coupling_summary <- group_summary(dat, "coupling")
  state_summary <- group_summary(dat, "drought_state")
  response <- fit_gpp_response(dat, bootstrap_n, block_n)
  sensitivity <- window_layer_sensitivity(response$data, cor_cut)
  notes <- drought_analysis_notes(met, gpp, matched, threshold, window_n,
                                  cor_cut, bootstrap_n, response)

  stopifnot(nrow(dat) == nrow(met), matched <= nrow(gpp),
            all(abs(met$r_SM1[is.finite(met$r_SM1)]) <= 1 + 1e-12),
            all(events$duration_days > 0),
            all(events$severity_z_days >= events$duration_days * threshold - 1e-10),
            all(response$coefficients$bootstrap_valid >= bootstrap_n * 0.8))

  write_drought_outputs(analysis_data, response, notes, threshold, cor_cut, out_dir)
  if (!is.null(plotter)) {
    plotter(dat = response$data, slopes = response$slopes,
            threshold = threshold, cor_cut = cor_cut,
            two_groups = response$two_groups, out_dir = file.path(out_dir, "images"))
  }

  cat("运行完成！结果目录：", out_dir, "\n")
  print(coupling_summary)
  print(state_summary)
  print(response$slopes)
  invisible(list(data = response$data, events = events,
                 coupling_summary = coupling_summary, state_summary = state_summary,
                 model = response$model, coefficients = response$coefficients,
                 mean_difference = response$mean_difference, slopes = response$slopes,
                 sensitivity = sensitivity, quality = inputs$quality, notes = notes))
}
