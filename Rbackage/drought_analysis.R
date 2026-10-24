# 旧入口，调用完整分析脚本
local({
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
  source(file.path(project_root, "scripts", "run_drought_analysis.R"),
         local = parent.env(environment()), encoding = "UTF-8")
})
