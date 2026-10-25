
seasonal_z <- function(x) {
  z <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) >= 5L && sd(x[ok]) > 0) {
    z[ok] <- (x[ok] - mean(x[ok])) / sd(x[ok])
  }
  z
}

event_table <- function(dates, index, days, label, threshold) {
  dry <- is.finite(index) & index >= threshold
  gap <- c(TRUE, as.numeric(diff(dates)) > 8)
  starts <- which(dry & (gap | !c(FALSE, head(dry, -1))))
  empty <- data.frame(index = character(), start = as.Date(character()),
                      end = as.Date(character()), n_periods = integer(),
                      duration_days = numeric(), intensity = numeric(),
                      severity_z_days = numeric(), boundary_or_gap = logical())
  if (!length(starts)) return(empty)
  events <- empty
  for (i in starts) {
    j <- i
    while (j < length(dry) && dry[j + 1L] && !gap[j + 1L]) j <- j + 1L
    ids <- i:j
    censored <- i == 1L || j == length(dry)
    if (i > 1L) censored <- censored || gap[i] || !is.finite(index[i - 1L])
    if (j < length(dry)) censored <- censored || gap[j + 1L] || !is.finite(index[j + 1L])
    event <- data.frame(index = label, start = dates[i], end = dates[j] + days[j] - 1,
               n_periods = length(ids), duration_days = sum(days[ids]),
               intensity = weighted.mean(index[ids], days[ids]),
               severity_z_days = sum(index[ids] * days[ids]),
               boundary_or_gap = censored)
    events <- rbind(events, event)
  }
  events
}

rolling_cor <- function(x, y, dates, width) {
  ans <- rep(NA_real_, length(x))
  if (length(x) < width) return(ans)
  for (i in seq.int(width, length(x))) {
    ids <- seq.int(i - width + 1L, i)
    if (any(as.numeric(diff(dates[ids])) > 8)) next
    ok <- is.finite(x[ids]) & is.finite(y[ids])
    if (sum(ok) >= ceiling(width * 0.8) &&
        sd(x[ids][ok]) > 0 && sd(y[ids][ok]) > 0) {
      ans[i] <- cor(x[ids][ok], y[ids][ok])
    }
  }
  ans
}

prepare_drought_indices <- function(met, threshold, window_n, cor_cut) {
  met |>
    mutate(slot = (as.integer(format(date, "%j")) - 1L) %/% 8L + 1L) |>
    group_by(slot) |>
    mutate(Z_VPD = seasonal_z(VPD), Z_SM1 = seasonal_z(SM1)) |>
    ungroup() |>
    arrange(date) |>
    mutate(
      D_VPD = Z_VPD, D_SM1 = -Z_SM1,
      days = c(pmin(8, as.numeric(diff(date))), 8),
      r_SM1 = rolling_cor(D_VPD, D_SM1, date, window_n),
      coupling = case_when(
        is.na(r_SM1) ~ NA_character_,
        r_SM1 >= cor_cut ~ "Synchronous",
        TRUE ~ "Asynchronous"),
      drought_state = case_when(
        !is.finite(D_VPD) | !is.finite(D_SM1) ~ NA_character_,
        D_VPD >= threshold & D_SM1 >= threshold ~ "Both",
        D_VPD >= threshold ~ "VPD_only",
        D_SM1 >= threshold ~ "SM_only",
        TRUE ~ "Neither")) |>
    as.data.frame()
}
