# 计算干旱指标和耦合类型

# 按年内8天时段，用全期数据计算季节性z分数
# SM指数取负号，数值越大越干；不是SPI/SPEI
seasonal_z <- function(x, slot) {
  z <- rep(NA_real_, length(x))
  for (s in unique(slot)) {
    ids <- which(slot == s & is.finite(x))
    if (length(ids) >= 5L && sd(x[ids]) > 0) {
      z[ids] <- (x[ids] - mean(x[ids])) / sd(x[ids])
    }
  }
  z
}

# 统计干旱强度、时长和烈度
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

# 向后滑动相关：至少80%有效配对，不跨越8天以上缺口
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
  met$slot <- (as.integer(format(met$date, "%j")) - 1L) %/% 8L + 1L

  for (v in c("VPD", "SM1")) met[[paste0("Z_", v)]] <- seasonal_z(met[[v]], met$slot)
  met$D_VPD <- met$Z_VPD
  met$D_SM1 <- -met$Z_SM1
  # 两种指数越大越干，正相关表示同向变化
  met$days <- c(pmin(8, as.numeric(diff(met$date))), 8)

  met$r_SM1 <- rolling_cor(met$D_VPD, met$D_SM1, met$date, window_n)
  # 按相关系数分类；耦合同步不等于同时干旱
  met$coupling <- ifelse(is.na(met$r_SM1), NA_character_,
                         ifelse(met$r_SM1 >= cor_cut, "Synchronous", "Asynchronous"))
  # 判断两种干旱是否同时发生
  met$drought_state <- NA_character_
  ok <- is.finite(met$D_VPD) & is.finite(met$D_SM1)
  a <- met$D_VPD[ok] >= threshold
  b <- met$D_SM1[ok] >= threshold
  met$drought_state[ok] <- ifelse(a & b, "Both", ifelse(a, "VPD_only",
                                                     ifelse(b, "SM_only", "Neither")))
  met
}

