# 干旱分析运行入口：支持Rscript和RStudio的Source。
# 终端（项目根目录）：Rscript scripts/run_drought_analysis.R
# R控制台（项目根目录）：source("scripts/run_drought_analysis.R", encoding = "UTF-8")
# 完整结果保存为drought_result；数据和输出路径根据本文件位置确定。
drought_result <- local({
  # source()优先使用当前文件的ofile；直接Rscript执行时读取--file。
  # chdir=TRUE时相对ofile可能已失效，因此同时检查当前目录下的文件名。
  source_files <- lapply(sys.frames(), function(frame) frame$ofile)
  source_files <- Filter(function(path) is.character(path) && length(path) == 1L,
                         source_files)
  if (length(source_files)) {
    entry_file <- tail(source_files, 1L)[[1L]]
  } else {
    file_arg <- commandArgs(trailingOnly = FALSE)
    file_arg <- file_arg[startsWith(file_arg, "--file=")]
    if (!length(file_arg)) stop("请通过Rscript或source()运行本入口文件")
    entry_file <- sub("^--file=", "", file_arg[1L])
  }
  if (!file.exists(entry_file) && file.exists(basename(entry_file))) {
    entry_file <- basename(entry_file)
  }
  entry_file <- normalizePath(entry_file, mustWork = TRUE)
  project_root <- dirname(dirname(entry_file))
  sys.source(file.path(project_root, "src", "drought_analysis.R"),
             envir = environment())
  sys.source(file.path(project_root, "scripts", "codeforimages", "plot_drought_results.R"),
             envir = environment())

  # 所有运行参数集中在这里；修改阈值时需保留其方法学含义。
  config <- list(
    data_dir = file.path(project_root, "Rbackage"),
    out_dir = file.path(project_root, "output"),
    threshold = 1,          # 干旱阈值：季节性z分数>=1
    window_n = 23L,        # 向后23期滑动窗口，约半年
    cor_cut = 0.3,         # 同步耦合的操作性相关阈值
    bootstrap_n = 999L,
    block_n = 23L,
    seed = 20261023,
    plotter = plot_drought_results # 图片统一保存到output/images/
  )
  do.call(run_drought_analysis, config)
})
