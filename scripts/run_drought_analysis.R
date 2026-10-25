drought_result <- local({
  suppressPackageStartupMessages(library(dplyr))
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
  modules <- c("data_io.R", "drought_indices.R", "gpp_response.R", "drought_analysis.R")
  for (module in modules) {
    sys.source(file.path(project_root, "src", module), envir = environment())
  }
  sys.source(file.path(project_root, "scripts", "codeforimages", "plot_drought_results.R"),
             envir = environment())

  data_dir <- file.path(project_root, "Rbackage")
  out_dir <- file.path(project_root, "output")
  threshold <- 1
  window_n <- 23L
  cor_cut <- 0.3
  result <- run_drought_analysis(data_dir, out_dir, threshold, window_n, cor_cut)
  plot_drought_results(result$data, file.path(out_dir, "images"))
  result
})
