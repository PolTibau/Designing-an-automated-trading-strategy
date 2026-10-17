## Full analysis: back-tests, rolling windows, robustness checks, tables and figures.
## Usage (from anywhere):  Rscript R/run_analysis.R
## Needs data/DB.csv and data/GOOG.csv (see README for the format). Writes to results/.

library(ggplot2)

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- if (length(script)) normalizePath(file.path(dirname(script), "..")) else getwd()
source(file.path(root, "R", "strategies.R"))

out_dir <- file.path(root, "results")
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------- parameters ----
## Chosen once, a priori (taken from the original script); NOT tuned on results.
lim <- 10          # EMA length for the news-flow filter
longshort <- -1    # -1: allow short positions; 0: go flat instead of short
thr <- 0.005       # |Trend| threshold for strategy 4 (0.5% of the open)
tcost <- 0.001     # proportional cost per unit traded (10 bps)
tickers <- c(DB = "Deutsche Bank", GOOG = "Google")
labels <- c(S1 = "1: Sentiment (news spikes)", S2 = "2: Sentiment + PP (news > EMA)",
            S3 = "3: Sentiment + PP", S4 = "4: Pivot points only", BH = "Buy & Hold")
cols <- c("1: Sentiment (news spikes)" = "#e69f00", "2: Sentiment + PP (news > EMA)" = "#009e73",
          "3: Sentiment + PP" = "#0072b2", "4: Pivot points only" = "#d55e00",
          "Buy & Hold" = "grey35")

## Positions held during each day for the four strategies + buy & hold.
build_exposures <- function(f, size = 1, trend = f$Trend) {
  list(
    S1 = lag_position(signal_1(f, lim, longshort, Pl = 65, Nl = 35)),
    S2 = lag_position(signal_2(f, lim, longshort, Pl = 65, Nl = 35)),
    S3 = lag_position(signal_3(f, lim, longshort, Pl = 50, Nl = 50)),
    S4 = lag_position(signal_4(f, longshort, thr, size, trend)),
    BH = rep(1, nrow(f))
  )
}

run_all <- function(f, exposures, cost) {
  lapply(exposures, function(e) backtest(f$AdjClose, e, cost))
}

metrics_table <- function(bts) {
  m <- sapply(names(bts), function(k) performance(bts[[k]]$ret, bts[[k]]$n_trades))
  colnames(m) <- labels[colnames(m)]
  m
}

## ---------------------------------------------------------------- run ----
curves <- list(); dds <- list(); full <- list(); roll <- list()
leak <- list(); cost_sens <- list(); halves <- list(); lever <- list()

for (tk in names(tickers)) {
  f <- build_features(load_market_data(file.path(root, "data", paste0(tk, ".csv"))))
  expo <- build_exposures(f)
  bts <- run_all(f, expo, tcost)

  ## Full-period metrics
  m <- metrics_table(bts)
  full[[tk]] <- data.frame(Ticker = tk, Metric = rownames(m), m, check.names = FALSE)

  ## Equity curves and drawdowns
  eq <- do.call(rbind, lapply(names(bts), function(k) {
    e <- cumprod(1 + bts[[k]]$ret)
    data.frame(Ticker = tickers[[tk]], Date = f$Date, Strategy = labels[[k]],
               Equity = e, Drawdown = e / cummax(e) - 1)
  }))
  curves[[tk]] <- eq

  ## Rolling windows (126 and 252 days)
  for (w in c(126, 252)) {
    per <- lapply(names(bts), function(k) {
      p <- rolling_performance(bts[[k]]$ret, w)
      data.frame(Ticker = tickers[[tk]], Window = w, Strategy = labels[[k]], p,
                 check.names = FALSE)
    })
    roll[[paste(tk, w)]] <- do.call(rbind, per)
  }

  ## Look-ahead bias: strategy 4 with the original (next-day) trend vs the corrected one
  leaky <- backtest(f$AdjClose, lag_position(signal_4(f, longshort, thr, 1, leaky_trend(f))), tcost)
  fixed <- bts$S4
  leak[[tk]] <- rbind(
    data.frame(Ticker = tickers[[tk]], Date = f$Date, Version = "Original (uses next-day pivot point)",
               Equity = cumprod(1 + leaky$ret)),
    data.frame(Ticker = tickers[[tk]], Date = f$Date, Version = "Corrected (uses information up to t)",
               Equity = cumprod(1 + fixed$ret))
  )

  ## Transaction-cost sensitivity (final cumulative return)
  costs <- c(0, 0.0005, 0.001, 0.0025, 0.005)
  cost_sens[[tk]] <- do.call(rbind, lapply(costs, function(cst) {
    b <- run_all(f, expo, cst)
    data.frame(Ticker = tickers[[tk]], tcost = cst, Strategy = unname(labels[names(b)]),
               CumReturn = sapply(b, function(z) prod(1 + z$ret) - 1))
  }))

  ## First half vs second half of the sample (fixed parameters)
  mid <- floor(nrow(f) / 2)
  idx <- list("First half" = seq_len(mid), "Second half" = (mid + 1):nrow(f))
  halves[[tk]] <- do.call(rbind, lapply(names(idx), function(h) {
    do.call(rbind, lapply(names(bts), function(k) {
      p <- performance(bts[[k]]$ret[idx[[h]]])
      data.frame(Ticker = tk, Period = h,
                 From = format(min(f$Date[idx[[h]]])), To = format(max(f$Date[idx[[h]]])),
                 Strategy = labels[[k]], CumReturn = p[["Cumulative Return"]],
                 Sharpe = p[["Annualized Sharpe Ratio"]], MaxDD = p[["Maximum Drawdown"]])
    }))
  }))

  ## Leverage of strategy 4 (the "aggressiveness" knob)
  lever[[tk]] <- do.call(rbind, lapply(c(1, 2, 3, 5), function(s) {
    b <- backtest(f$AdjClose, lag_position(signal_4(f, longshort, thr, s)), tcost)
    p <- performance(b$ret, b$n_trades)
    data.frame(Ticker = tk, Leverage = s, CumReturn = p[["Cumulative Return"]],
               Sharpe = p[["Annualized Sharpe Ratio"]], Volatility = p[["Annualized Volatility"]],
               MaxDD = p[["Maximum Drawdown"]])
  }))
}

bind <- function(l) { x <- do.call(rbind, l); rownames(x) <- NULL; x }
curves <- bind(curves); leak <- bind(leak); cost_sens <- bind(cost_sens)
roll <- bind(roll); full <- bind(full); halves <- bind(halves); lever <- bind(lever)

## ------------------------------------------------------------ tables ----
write.csv(full, file.path(out_dir, "metrics_full_period.csv"), row.names = FALSE)
write.csv(halves, file.path(out_dir, "metrics_by_half.csv"), row.names = FALSE)
write.csv(lever, file.path(out_dir, "strategy4_leverage.csv"), row.names = FALSE)
write.csv(cost_sens, file.path(out_dir, "cost_sensitivity.csv"), row.names = FALSE)
roll_mean <- aggregate(. ~ Ticker + Window + Strategy, data = roll, FUN = mean,
                       na.rm = TRUE, na.action = na.pass)
n_win <- aggregate(`Annual Return` ~ Ticker + Window + Strategy, data = roll, FUN = length)
names(n_win)[4] <- "NumWindows"
roll_mean <- merge(roll_mean, n_win, by = c("Ticker", "Window", "Strategy"))
write.csv(roll_mean, file.path(out_dir, "metrics_rolling_mean.csv"), row.names = FALSE)

## ----------------------------------------------------------- figures ----
theme_set(theme_minimal(base_size = 11) +
            theme(legend.position = "bottom", legend.title = element_blank(),
                  plot.title.position = "plot"))
save <- function(p, name, w = 9, h = 5) {
  ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 150, bg = "white")
}
curves$Strategy <- factor(curves$Strategy, levels = labels)

save(ggplot(curves, aes(Date, Equity, colour = Strategy)) + geom_line(linewidth = 0.6) +
       facet_wrap(~Ticker, ncol = 1, scales = "free_y") + scale_colour_manual(values = cols) +
       guides(colour = guide_legend(nrow = 2)) +
       labs(title = "Growth of 1 unit invested (net of costs)", y = "Equity", x = NULL),
     "01_equity_curves.png", h = 7)

save(ggplot(curves, aes(Date, Drawdown, colour = Strategy)) + geom_line(linewidth = 0.5) +
       facet_wrap(~Ticker, ncol = 1) + scale_colour_manual(values = cols) +
       scale_y_continuous(labels = function(v) paste0(round(100 * v), "%")) +
       guides(colour = guide_legend(nrow = 2)) +
       labs(title = "Drawdowns", y = "Drawdown from previous peak", x = NULL),
     "02_drawdowns.png", h = 7)

roll$Strategy <- factor(roll$Strategy, levels = labels)
save(ggplot(subset(roll, Window == 252), aes(Strategy, `Annual Return`, fill = Strategy)) +
       geom_hline(yintercept = 0, colour = "grey50") + geom_boxplot(outlier.size = 0.6) +
       facet_wrap(~Ticker) + scale_fill_manual(values = cols) +
       scale_y_continuous(labels = function(v) paste0(round(100 * v), "%")) +
       theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
       guides(fill = guide_legend(nrow = 2)) +
       labs(title = "Annualised return across all 252-day rolling windows", y = "Annual return"),
     "03_rolling_252_returns.png", h = 5)

save(ggplot(subset(roll, Window == 126), aes(Strategy, `Annualized Sharpe Ratio`, fill = Strategy)) +
       geom_hline(yintercept = 0, colour = "grey50") + geom_boxplot(outlier.size = 0.6) +
       facet_wrap(~Ticker) + scale_fill_manual(values = cols) +
       theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
       guides(fill = guide_legend(nrow = 2)) +
       labs(title = "Sharpe ratio across all 126-day rolling windows", y = "Annualised Sharpe"),
     "04_rolling_126_sharpe.png", h = 5)

save(ggplot(leak, aes(Date, Equity, colour = Version)) + geom_line(linewidth = 0.7) +
       facet_wrap(~Ticker, ncol = 1, scales = "free_y") +
       scale_y_log10() + scale_colour_manual(values = c("#d55e00", "#0072b2")) +
       labs(title = "Strategy 4: effect of the look-ahead bug in the original script",
            y = "Equity (log scale)", x = NULL),
     "05_lookahead_bias.png", h = 6.5)

cost_sens$Strategy <- factor(cost_sens$Strategy, levels = labels)
save(ggplot(cost_sens, aes(tcost * 1e4, CumReturn, colour = Strategy)) +
       geom_line(linewidth = 0.7) + geom_point() + facet_wrap(~Ticker, scales = "free_y") +
       scale_colour_manual(values = cols) +
       scale_y_continuous(labels = function(v) paste0(round(100 * v), "%")) +
       guides(colour = guide_legend(nrow = 2)) +
       labs(title = "Sensitivity to transaction costs", x = "Cost per unit traded (bps)",
            y = "Cumulative return"),
     "06_cost_sensitivity.png", h = 5)

## ---------------------------------------------------------- console ----
options(width = 200, digits = 3)
cat("\n== Full period (net of", tcost * 1e4, "bps) ==\n"); print(full, row.names = FALSE)
cat("\n== Mean over rolling windows ==\n"); print(roll_mean, row.names = FALSE)
cat("\n== First vs second half ==\n"); print(halves, row.names = FALSE)
cat("\n== Strategy 4 leverage ==\n"); print(lever, row.names = FALSE)
cat("\nFigures and tables written to", out_dir, "\n")
