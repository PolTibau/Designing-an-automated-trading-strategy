## Trading signals, back-test engine and performance metrics.
## Pure functions only: no file reading, no plotting (see run_analysis.R).
##
## Timing convention (avoids look-ahead bias)
##   * Every signal is computed with information available at the CLOSE of day t.
##   * The resulting position is held during day t+1 (see lag_position()).
##   * Returns are close-to-close returns of the adjusted close.

library(TTR) # EMA

## ---------------------------------------------------------------- data ----

load_market_data <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  x <- x[, names(x) != "X", drop = FALSE] # unnamed index column written by the export
  x$Date <- as.Date(x$Date)
  x <- x[order(x$Date), ]
  rownames(x) <- NULL
  x
}

## Builds the indicators used by the strategies.
##   BBr   Bull/Bear ratio in [0, 100]; 50 = neutral (also when there is no news).
##   RVT   relative volume of text (news flow); 0 on days without news.
##   Trend (PivotPoint_t - Open_t) / Open_t, with PivotPoint_t = (High+Low+Close)/3.
##         Scale-free, so the same threshold works for a 10 EUR and a 1000 USD stock.
build_features <- function(x) {
  fill0 <- function(v) ifelse(is.na(v), 0, v)
  bull <- fill0(x$positivePartscr) + fill0(x$certaintyPartscr) + fill0(x$finupPartscr)
  bear <- fill0(x$negativePartscr) + fill0(x$uncertaintyPartscr) + fill0(x$findownPartscr)
  bbr <- 100 * bull / (bull + bear)
  bbr[is.na(bbr)] <- 50 # bull = bear = 0 -> no information
  pp <- (x$High + x$Low + x$Close) / 3
  data.frame(
    Date = x$Date,
    AdjClose = x$Adj.Close,
    BBr = bbr,
    RVT = fill0(x$RVT),
    Trend = (pp - x$Open) / x$Open,
    PP = pp,
    Open = x$Open
  )
}

## Trend as it was (incorrectly) built in the original script: it used the NEXT day's
## pivot point, which contains the next day's close. Kept ONLY to quantify the bias.
leaky_trend <- function(feat) {
  (c(feat$PP[-1], NA) - feat$Open) / feat$Open
}

## ------------------------------------------------------------- signals ----

## Position state machine: enter long / enter short (or go flat when longshort = 0)
## and HOLD the position until the opposite signal appears.
hold_position <- function(go_long, go_short, longshort = -1, size = 1) {
  go_long[is.na(go_long)] <- FALSE
  go_short[is.na(go_short)] <- FALSE
  pos <- numeric(length(go_long))
  cur <- 0
  for (i in seq_along(pos)) {
    if (go_long[i]) {
      cur <- size
    } else if (go_short[i]) {
      cur <- if (longshort == -1) -size else 0
    }
    pos[i] <- cur
  }
  pos
}

## Volume-of-news filter: RVT above `mult` times its EMA (needs more than `lim` days).
rvt_filter <- function(rvt, lim, mult = 1) {
  ok <- seq_along(rvt) > lim & rvt > mult * EMA(rvt, lim)
  ok[is.na(ok)] <- FALSE
  ok
}

## 1: sentiment only, when the news flow spikes (RVT > 1.5 * EMA).
signal_1 <- function(f, lim, longshort, Pl, Nl) {
  gate <- rvt_filter(f$RVT, lim, 1.5)
  hold_position(gate & f$BBr > Pl, gate & f$BBr < Nl, longshort)
}

## 2: sentiment + pivot-point trend, when news flow is above its EMA.
signal_2 <- function(f, lim, longshort, Pl, Nl) {
  gate <- rvt_filter(f$RVT, lim, 1)
  hold_position(gate & f$BBr > Pl & f$Trend > 0,
                gate & f$BBr < Nl & f$Trend < 0, longshort)
}

## 3: like 2 but without the news-flow filter (fewer restrictions to trade).
signal_3 <- function(f, lim, longshort, Pl, Nl) {
  hold_position(f$BBr > Pl & f$Trend > 0, f$BBr < Nl & f$Trend < 0, longshort)
}

## 4: pivot points only. `size` is the leverage: the aggressiveness knob.
signal_4 <- function(f, longshort, thr, size = 1, trend = f$Trend) {
  hold_position(trend > thr, trend < -thr, longshort, size)
}

## Position decided at the close of t is held during t+1.
lag_position <- function(position) c(0, position[-length(position)])

## ----------------------------------------------------------- back-test ----

## `exposure` is the position held DURING each day (already lagged).
## Proportional cost: `tcost` * |change in exposure|; the final position is liquidated.
backtest <- function(adj_close, exposure, tcost = 0) {
  n <- length(adj_close)
  ret <- c(0, adj_close[-1] / adj_close[-n] - 1)
  turnover <- abs(diff(c(0, exposure)))
  turnover[n] <- turnover[n] + abs(exposure[n]) # liquidation at the end
  net <- exposure * ret - tcost * turnover
  list(ret = net, n_trades = sum(turnover > 0))
}

## ------------------------------------------------------------ metrics ----

performance <- function(r, n_trades = NA, periods = 252) {
  equity <- cumprod(1 + r)
  cum <- equity[length(equity)] - 1
  ann_ret <- if (1 + cum > 0) (1 + cum)^(periods / length(r)) - 1 else -1
  ann_vol <- sd(r) * sqrt(periods)
  traded <- r[r != 0]
  dd <- equity / cummax(c(1, equity))[-1] - 1
  runs <- rle(dd < 0)
  c(
    "Cumulative Return" = cum,
    "Annual Return" = ann_ret,
    "Annualized Sharpe Ratio" = if (ann_vol > 0) ann_ret / ann_vol else NA_real_,
    "Win %" = if (length(traded)) mean(traded > 0) else NA_real_,
    "Annualized Volatility" = ann_vol,
    "Maximum Drawdown" = min(dd),
    "Max Length Drawdown" = if (any(runs$values)) max(runs$lengths[runs$values]) else 0,
    "n.trades" = n_trades
  )
}

## Performance in every window of `w` consecutive days (step = 1 day).
rolling_performance <- function(r, w) {
  n_windows <- length(r) - w + 1
  if (n_windows < 1) stop("Window size too large")
  metrics <- names(performance(r))[1:7]
  out <- t(vapply(seq_len(n_windows),
                  function(i) unname(performance(r[i:(i + w - 1)])[1:7]),
                  numeric(7)))
  colnames(out) <- metrics
  out
}
