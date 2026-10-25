local({
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
  dat <- read.csv(file.path(out_dir, "04_analysis_data.csv"), stringsAsFactors = FALSE,
                  na.strings = c("", "NA"))
  dat$date <- as.Date(dat$date)
  image_dir <- file.path(out_dir, "images")
  plot_drought_results(dat, image_dir)
  cat("绘图完成！图片目录：", image_dir, "\n")
})
