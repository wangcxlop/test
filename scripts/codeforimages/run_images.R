# 独立绘图入口：仅读取已有CSV和绘图参数，不重新计算分析或bootstrap。
# 项目根目录运行：Rscript scripts/codeforimages/run_images.R
# RStudio：source("scripts/codeforimages/run_images.R", encoding = "UTF-8")
local({
  # 根据入口位置定位项目，支持不同工作目录和source(chdir=TRUE)。
  source_files <- lapply(sys.frames(), function(frame) frame$ofile)
  source_files <- Filter(function(path) is.character(path) && length(path) == 1L,
                         source_files)
  if (length(source_files)) {
    entry_file <- tail(source_files, 1L)[[1L]]
  } else {
    file_arg <- commandArgs(trailingOnly = FALSE)
    file_arg <- file_arg[startsWith(file_arg, "--file=")]
    if (!length(file_arg)) stop("请通过Rscript或source()运行本绘图入口")
    entry_file <- sub("^--file=", "", file_arg[1L])
  }
  if (!file.exists(entry_file) && file.exists(basename(entry_file))) {
    entry_file <- basename(entry_file)
  }
  entry_file <- normalizePath(entry_file, mustWork = TRUE)
  project_root <- dirname(dirname(dirname(entry_file)))
  sys.source(file.path(dirname(entry_file), "plot_drought_results.R"),
             envir = environment())
  out_dir <- file.path(project_root, "output")
  inputs <- read_plot_inputs(out_dir)
  image_dir <- file.path(out_dir, "images")
  do.call(plot_drought_results, c(inputs, list(out_dir = image_dir)))
  cat("绘图完成！图片目录：", image_dir, "\n")
})
