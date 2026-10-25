
paired_correlation <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  r <- if (n >= 3L && sd(x[ok]) > 0 && sd(y[ok]) > 0) cor(x[ok], y[ok]) else NA_real_
  list(n = n, r = r)
}

group_summary <- function(dat, col, groups) {
  dat |>
    rename(group = all_of(col)) |>
    mutate(group = factor(group, levels = groups)) |>
    filter(!is.na(group), is.finite(GPP)) |>
    group_by(group, .drop = FALSE) |>
    summarise(
      n_GPP = n(),
      mean_GPP = if (n_GPP > 0) mean(GPP) else NA_real_,
      n_VPD = sum(is.finite(D_VPD)),
      r_VPD = paired_correlation(GPP, D_VPD)$r,
      n_SM1 = sum(is.finite(D_SM1)),
      r_SM1 = paired_correlation(GPP, D_SM1)$r,
      .groups = "drop") |>
    mutate(grouping = col, group = as.character(group)) |>
    relocate(grouping, group) |>
    as.data.frame()
}

analyze_gpp_response <- function(dat) {
  bind_rows(
    group_summary(dat, "coupling", c("Asynchronous", "Synchronous")),
    group_summary(dat, "drought_state", c("Both", "VPD_only", "SM_only", "Neither")))
}
