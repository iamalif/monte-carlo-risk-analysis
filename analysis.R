# Risk Analysis AMI22Z, Home Assignment 2, Group 15
# Run with Rscript analysis.R or RStudio's Source button. The script reads the frozen Yahoo files in raw/,
# writes result tables beside itself, and saves the plots in figures/.
# The price files were downloaded on 21 October 2026.

# Look beside the script first, then in the working directory. This covers both
# Rscript and RStudio's Source button without asking anyone to edit paths.
args <- commandArgs(trailingOnly = FALSE)
script_files <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
source_files <- as.character(unlist(lapply(sys.frames(), function(frame) {
  if (is.character(frame$ofile)) frame$ofile else character(0)
}), use.names = FALSE))
candidate_dirs <- unique(normalizePath(
  c(dirname(script_files), dirname(rev(source_files)), getwd()),
  mustWork = FALSE
))
has_inputs <- vapply(candidate_dirs, function(path) {
  all(file.exists(file.path(path, "raw", c("INTC.csv", "SSAB-A.ST.csv", "SEK_X.csv"))))
}, logical(1))
if (!any(has_inputs)) {
  stop("Cannot find the raw/ CSV files. Extract the full ZIP, then run analysis.R from that folder.")
}
script_dir <- candidate_dirs[which(has_inputs)[1L]]
raw_dir <- file.path(script_dir, "raw")
fig_dir <- file.path(script_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Yahoo has a few empty FX rows. They are dropped before we look up dates.
read_prices <- function(symbol) {
  x <- read.csv(file.path(raw_dir, paste0(gsub("=", "_", symbol), ".csv")))
  x$date <- as.Date(x$date)
  x <- x[order(x$date), ]
  x <- x[!is.na(x$close) & x$close > 0, ]
  stopifnot(!anyDuplicated(x$date))
  x
}

# Consecutive rows mean consecutive trading dates for that share, not calendar days.
returns <- function(x) {
  data.frame(date = x$date[-1], r = diff(log(x$close)))
}

# Wilson intervals behave better than the usual normal interval for proportions.
wilson <- function(success, n, level = 0.95) {
  z <- qnorm((1 + level) / 2)
  p <- success / n
  den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(lower = max(0, center - half), upper = min(1, center + half))
}

mean_ci <- function(x) {
  h <- qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
  c(lower = mean(x) - h, upper = mean(x) + h)
}

# This chi-square interval is exact only under independent normal returns.
sd_ci <- function(x) {
  df <- length(x) - 1L
  ss <- df * var(x)
  sqrt(c(lower = ss / qchisq(0.975, df),
         upper = ss / qchisq(0.025, df)))
}

# This extra check is useful because the historical returns have heavy tails.
# Resampling 20-day stretches keeps some of the volatility clustering; wrapping
# at the end lets any trading date start a block. The interval still assumes
# the pre-2022 return process is reasonably stable over the sample.
excess_kurtosis <- function(x) {
  centered <- x - mean(x)
  mean(centered^4) / mean(centered^2)^2 - 3
}

block_sd_ci <- function(x, seed, replicates = 10000L, block_length = 20L) {
  set.seed(seed)
  n <- length(x)
  blocks_needed <- ceiling(n / block_length)
  estimates <- numeric(replicates)
  for (b in seq_len(replicates)) {
    starts <- sample.int(n, blocks_needed, replace = TRUE)
    indices <- as.vector(t(outer(starts, 0:(block_length - 1L),
                                 function(a, offset) ((a + offset - 1L) %% n) + 1L)))
    estimates[b] <- sd(x[indices[seq_len(n)]])
  }
  unname(quantile(estimates, c(0.025, 0.975)))
}

# Problem 1. Only returns available before the investment go into the fit.
intel <- read_prices("INTC")
ssab <- read_prices("SSAB-A.ST")
fx <- read_prices("SEK=X") # SEK per USD; used only for purchase/cash-out conversion.
cutoff <- as.Date("2021-12-31")
intel_r <- returns(intel)
ssab_r <- returns(ssab)
intel_train <- intel_r[intel_r$date <= cutoff, ]
ssab_train <- ssab_r[ssab_r$date <= cutoff, ]
paired <- merge(intel_train, ssab_train, by = "date", suffixes = c("_intc", "_ssab"))
# Each share uses its full history for mean and volatility. Correlation uses
# only dates when both markets have a return.
mu <- c(INTC = mean(intel_train$r), SSAB_A = mean(ssab_train$r))
sigma <- c(INTC = sd(intel_train$r), SSAB_A = sd(ssab_train$r))
rho <- cor(paired$r_intc, paired$r_ssab)
cov_daily <- matrix(c(sigma[1]^2, rho * sigma[1] * sigma[2],
                      rho * sigma[1] * sigma[2], sigma[2]^2), 2, 2)
initial_date <- as.Date("2022-01-03")
end_date <- as.Date("2026-06-30")
price_on <- function(x, date) {
  value <- x$close[x$date == date]
  if (length(value) != 1L) stop("Missing date ", date)
  value
}
adjusted_on <- function(x, date) {
  value <- x$adjusted_close[x$date == date]
  if (length(value) != 1L || is.na(value)) stop("Missing adjusted close on ", date)
  value
}
p0 <- c(INTC = price_on(intel, initial_date), SSAB_A = price_on(ssab, initial_date))
p_actual <- c(INTC = price_on(intel, end_date), SSAB_A = price_on(ssab, end_date))
adj0 <- c(INTC = adjusted_on(intel, initial_date),
          SSAB_A = adjusted_on(ssab, initial_date))
adj_actual <- c(INTC = adjusted_on(intel, end_date),
                SSAB_A = adjusted_on(ssab, end_date))
fx0 <- price_on(fx, initial_date)
fx_end <- price_on(fx, end_date)
horizon <- as.integer(4.5 * 252)
n_paths <- 20000L # At least 10,000 full daily trajectories, per instructor guidance.
batch <- 500L
# Batches keep the full-path simulation from holding every daily price in memory.
terminal <- numeric(n_paths)
terminal_intel <- numeric(n_paths)
terminal_ssab <- numeric(n_paths)
sample_paths <- NULL
sample_intel <- NULL
sample_ssab <- NULL
set.seed(150103)
root <- chol(cov_daily)
for (first in seq.int(1L, n_paths, by = batch)) {
  last <- min(n_paths, first + batch - 1L)
  k <- last - first + 1L
  # The Cholesky factor gives the two shares their fitted correlation.
  innovation <- matrix(rnorm(horizon * k * 2L), ncol = 2L) %*% root
  ri <- matrix(innovation[, 1] + mu[1], nrow = horizon, ncol = k)
  rs <- matrix(innovation[, 2] + mu[2], nrow = horizon, ncol = k)
  end_i <- p0[1] * exp(colSums(ri))
  end_s <- p0[2] * exp(colSums(rs))
  terminal_intel[first:last] <- end_i
  terminal_ssab[first:last] <- end_s
  terminal[first:last] <- 500 * end_i / p0[1] + 500 * end_s / p0[2]
  # Only the first paths are kept in memory for the figures. All 20,000 paths
  # are counted in the portfolio results.
  if (is.null(sample_paths)) {
    show <- 1:60
    sample_intel <- p0[1] * exp(apply(ri[, show, drop = FALSE], 2, cumsum))
    sample_ssab <- p0[2] * exp(apply(rs[, show, drop = FALSE], 2, cumsum))
    sample_paths <- cbind(
      500 * sample_intel / p0[1] + 500 * sample_ssab / p0[2]
    )
  }
}
profit <- terminal - 1000
p_profitable <- mean(terminal > 1000)
profit_ci <- mean_ci(profit)
p_profitable_ci <- wilson(sum(terminal > 1000), n_paths)
# The main comparison keeps the purchase FX rate fixed. The cash-out result
# below converts the SSAB holding back to dollars at the final observed rate.
actual_constant_fx <- 500 * p_actual[1] / p0[1] + 500 * p_actual[2] / p0[2]
actual_usd_cashout <- 500 * p_actual[1] / p0[1] +
  500 * p_actual[2] / p0[2] * fx0 / fx_end
# Adj Close ratios are a reinvestment proxy, not a ledger of cash dividends.
adjusted_components <- 500 * adj_actual / adj0
actual_adjusted_constant_fx <- sum(adjusted_components)
actual_adjusted_usd_cashout <- adjusted_components[1] +
  adjusted_components[2] * fx0 / fx_end
actual_percentile <- mean(terminal <= actual_constant_fx)
analytic_mean <- sum(500 * exp(horizon * (mu + sigma^2 / 2)))
drift_factor_ci <- c(
  exp(horizon * mean_ci(intel_train$r)),
  exp(horizon * mean_ci(ssab_train$r))
)
bootstrap_sd_intel <- block_sd_ci(intel_train$r, seed = 3141L)
bootstrap_sd_ssab <- block_sd_ci(ssab_train$r, seed = 3142L)
purchase_shares <- c(INTC = 500 / p0[1], SSAB_A = 500 * fx0 / p0[2])
actual_return_steps <- c(INTC = sum(intel$date > initial_date & intel$date <= end_date),
                         SSAB_A = sum(ssab$date > initial_date & ssab$date <= end_date))
param_stats <- data.frame(
  asset = c("INTC", "SSAB-A.ST"),
  observations = c(nrow(intel_train), nrow(ssab_train)),
  mean_log_return = as.numeric(mu),
  mean_ci_lower = c(mean_ci(intel_train$r)[1], mean_ci(ssab_train$r)[1]),
  mean_ci_upper = c(mean_ci(intel_train$r)[2], mean_ci(ssab_train$r)[2]),
  sd_log_return = as.numeric(sigma),
  sd_ci_lower = c(sd_ci(intel_train$r)[1], sd_ci(ssab_train$r)[1]),
  sd_ci_upper = c(sd_ci(intel_train$r)[2], sd_ci(ssab_train$r)[2]),
  excess_kurtosis = c(excess_kurtosis(intel_train$r),
                      excess_kurtosis(ssab_train$r)),
  block_bootstrap_sd_ci_lower = c(bootstrap_sd_intel[1], bootstrap_sd_ssab[1]),
  block_bootstrap_sd_ci_upper = c(bootstrap_sd_intel[2], bootstrap_sd_ssab[2]),
  variance_log_return = as.numeric(sigma^2),
  initial_close = as.numeric(p0),
  final_close = as.numeric(p_actual),
  currency = c("USD", "SEK")
)
write.csv(param_stats, file.path(script_dir, "stock_parameters.csv"), row.names = FALSE)
stock_results <- data.frame(
  metric = c("number_of_paths", "trading_days_simulated", "paired_training_days",
             "return_correlation", "expected_final_value_constant_fx",
             "expected_profit_constant_fx", "expected_profit_ci_lower",
             "expected_profit_ci_upper", "probability_profitable_constant_fx",
             "probability_profitable_ci_lower", "probability_profitable_ci_upper",
             "median_final_value_constant_fx", "final_value_p05", "final_value_p95",
             "actual_final_value_constant_fx", "actual_profit_constant_fx",
             "actual_forecast_percentile", "purchase_usdsek", "final_usdsek",
             "actual_final_value_usd_cashout", "actual_profit_usd_cashout",
             "analytic_expected_final_value", "purchase_intel_shares",
             "purchase_ssab_shares", "purchase_ssab_sek",
             "actual_adjusted_final_value_constant_fx",
             "actual_adjusted_final_value_usd_cashout",
             "actual_adjusted_intel_component",
             "actual_adjusted_ssab_component_constant_fx",
             "actual_intel_return_steps", "actual_ssab_return_steps",
             "histogram_count_above_5000",
             "intc_drift_factor_ci_lower", "intc_drift_factor_ci_upper",
             "ssab_drift_factor_ci_lower", "ssab_drift_factor_ci_upper"),
  value = c(n_paths, horizon, nrow(paired), rho, mean(terminal), mean(profit),
            profit_ci, p_profitable, p_profitable_ci,
            median(terminal), quantile(terminal, 0.05), quantile(terminal, 0.95),
            actual_constant_fx, actual_constant_fx - 1000, actual_percentile,
            fx0, fx_end, actual_usd_cashout, actual_usd_cashout - 1000,
            analytic_mean, purchase_shares, 500 * fx0,
            actual_adjusted_constant_fx, actual_adjusted_usd_cashout,
            adjusted_components, actual_return_steps, sum(terminal >= 5000),
            drift_factor_ci)
)
write.csv(stock_results, file.path(script_dir, "stock_results.csv"), row.names = FALSE)

# For the shorter fitting windows we only need final values. Under this normal
# model, one draw of the summed return has the same terminal distribution as
# 1,134 daily draws. These draws are separate from the main 20,000 paths.
set.seed(160101)
sensitivity_z <- matrix(rnorm(200000L), ncol = 2L)
window_sensitivity <- function(start, label, independent = FALSE) {
  prices_i <- intel[intel$date >= as.Date(start) & intel$date <= cutoff, ]
  prices_s <- ssab[ssab$date >= as.Date(start) & ssab$date <= cutoff, ]
  ri <- returns(prices_i)
  rs <- returns(prices_s)
  aligned <- merge(ri, rs, by = "date", suffixes = c("_i", "_s"))
  m <- c(mean(ri$r), mean(rs$r))
  s <- c(sd(ri$r), sd(rs$r))
  correlation <- if (independent) 0 else cor(aligned$r_i, aligned$r_s)
  covariance <- matrix(c(s[1]^2, correlation * prod(s),
                         correlation * prod(s), s[2]^2), 2, 2)
  cumulative <- sqrt(horizon) * sensitivity_z %*% chol(covariance) +
    matrix(horizon * m, nrow(sensitivity_z), 2L, byrow = TRUE)
  value <- 500 * exp(cumulative[, 1]) + 500 * exp(cumulative[, 2])
  data.frame(window = label, start_date = start,
             intc_daily_mean = m[1], ssab_daily_mean = m[2],
             return_correlation = correlation,
             analytic_mean_final_value = sum(500 * exp(horizon * (m + s^2 / 2))),
             simulated_profit_probability = mean(value > 1000),
             probability_ci_lower = wilson(sum(value > 1000), length(value))[1],
             probability_ci_upper = wilson(sum(value > 1000), length(value))[2],
             simulations = length(value))
}
sensitivity <- do.call(rbind, list(
  window_sensitivity("2009-10-01", "2009-2021"),
  window_sensitivity("2012-01-01", "2012-2021"),
  window_sensitivity("2017-01-01", "2017-2021"),
  window_sensitivity("2019-01-01", "2019-2021"),
  window_sensitivity("2009-10-01", "2009-2021, rho=0", TRUE)
))
write.csv(sensitivity, file.path(script_dir, "stock_sensitivity.csv"),
          row.names = FALSE)

hist_breaks <- seq(0, 5000, by = 500)
hist_counts <- hist(terminal[terminal < 5000], breaks = hist_breaks,
                    plot = FALSE, right = FALSE)$counts
write.csv(data.frame(bin = c(paste0(hist_breaks[-length(hist_breaks)], "-",
                                    hist_breaks[-1] - 1), "5000+"),
                     count = c(hist_counts, sum(terminal >= 5000))),
          file.path(script_dir, "stock_histogram.csv"), row.names = FALSE)

observed_i <- intel[intel$date >= initial_date & intel$date <= end_date, ]
observed_s <- ssab[ssab$date >= initial_date & ssab$date <= end_date, ]
observed_pair <- merge(observed_i[, c("date", "close")],
                       observed_s[, c("date", "close")], by = "date",
                       suffixes = c("_i", "_s"))
# The markets have different holidays. Plot the observed portfolio only where
# both closes exist, and scale calendar time to the simulated 1,134-step axis.
observed_pair$day <- as.numeric(observed_pair$date - initial_date) /
  as.numeric(end_date - initial_date) * horizon
observed_pair$value <- 500 * observed_pair$close_i / p0[1] +
  500 * observed_pair$close_s / p0[2]

png(file.path(fig_dir, "stock_history.png"), width = 1500, height = 850, res = 150)
par(mfrow = c(2, 1), mar = c(3.5, 5, 2, 1), oma = c(0, 0, 1, 0))
for (item in list(list(x = intel, name = "Intel closing price", unit = "USD"),
                  list(x = ssab, name = "SSAB-A closing price", unit = "SEK"))) {
  x <- item$x[item$x$date <= end_date, ]
  plot(x$date, x$close, type = "l", lwd = 1.2, col = "#245A84",
       xlab = "", ylab = paste0(item$unit, " per share"), main = item$name)
  abline(v = initial_date, col = "#B33C3C", lty = 2)
}
dev.off()

png(file.path(fig_dir, "stock_share_paths.png"), width = 1500, height = 1100, res = 150)
par(mfrow = c(2, 1), mar = c(4, 5, 2, 1))
for (asset in 1:2) {
  p <- p0[asset]
  m <- mu[asset]
  s <- sigma[asset]
  days <- 0:horizon
  lower <- p * exp(days * m + qnorm(0.05) * s * sqrt(days))
  middle <- p * exp(days * m)
  upper <- p * exp(days * m + qnorm(0.95) * s * sqrt(days))
  observed <- if (asset == 1) observed_i else observed_s
  paths <- if (asset == 1) sample_intel else sample_ssab
  names <- c("Intel (USD)", "SSAB-A (SEK)")
  plot(days, middle, type = "n", ylim = c(0, max(upper, observed$close) * 1.04),
       xlab = "Simulated trading day; actual dates scaled to horizon",
       ylab = names[asset], main = names[asset])
  polygon(c(days, rev(days)), c(lower, rev(upper)),
          col = adjustcolor("#AFC5DA", alpha.f = 0.55), border = NA)
  matlines(1:horizon, paths[, 1:15, drop = FALSE],
           col = adjustcolor("#6B83A0", alpha.f = 0.18), lty = 1)
  lines(days, middle, col = "#174263", lwd = 2)
  observed_days <- as.numeric(observed$date - initial_date) /
    as.numeric(end_date - initial_date) * horizon
  lines(observed_days, observed$close, col = "#B33C3C", lwd = 2)
  if (asset == 1) legend("topleft", c("5th-95th model band", "Model median",
                                      "Observed close"),
                         fill = c("#AFC5DA", NA, NA),
                         border = c(NA, NA, NA),
                         lty = c(NA, 1, 1),
                         col = c(NA, "#174263", "#B33C3C"), bty = "n", cex = 0.8)
}
dev.off()

png(file.path(fig_dir, "stock_paths.png"), width = 1500, height = 900, res = 150)
matplot(0:horizon, rbind(rep(1000, ncol(sample_paths)), sample_paths),
        type = "l", lty = 1, col = adjustcolor("#6B83A0", alpha.f = 0.28),
        xlab = "Simulated trading day", ylab = "Portfolio value (USD, fixed purchase FX)",
        main = "First 60 simulated portfolio trajectories and observed path")
abline(h = 1000, lty = 2, col = "#333333")
lines(observed_pair$day, observed_pair$value, col = "#B33C3C", lwd = 2.5)
legend("topleft", "Observed price-only path", col = "#B33C3C",
       lty = 1, lwd = 2.5, bty = "n")
dev.off()

png(file.path(fig_dir, "stock_distribution.png"), width = 1500, height = 900, res = 150)
hist(terminal, breaks = 90, xlim = c(0, quantile(terminal, 0.99)),
     col = "#AFC5DA", border = "white", xlab = "Final portfolio value (USD, fixed purchase FX)",
     main = "Simulated final value and observed price-only outcome")
abline(v = 1000, col = "#444444", lwd = 2, lty = 2)
abline(v = actual_constant_fx, col = "#B33C3C", lwd = 3)
legend("topright", legend = c("Initial $1,000", "Observed, fixed FX"),
       col = c("#444444", "#B33C3C"), lty = c(2, 1), lwd = c(2, 3), bty = "n")
dev.off()

# Problem 3. Discharges happen before the day's arrivals. A newly admitted
# patient remains for that day, which gives the geometric stay its one-day
# minimum. Anyone who cannot be admitted leaves rather than waiting.
arrival_rate <- 14.8
leave_probability <- 1 / 2.9
survive_probability <- 1 - leave_probability
burn_days <- 365L
study_days <- 365L
n_hospital_paths <- 10000L

simulate_hospital <- function(capacity, n_paths = n_hospital_paths,
                              burn = burn_days, days = study_days, seed = 1000 + capacity) {
  set.seed(seed)
  occupied <- integer(n_paths)
  lost_total <- numeric(n_paths)
  arrival_total <- numeric(n_paths)
  within_total <- numeric(n_paths)
  occupied_total <- numeric(n_paths)
  # The 85% rule is about an integer number of occupied beds.
  limit <- floor(0.85 * capacity + 1e-10)
  for (day in seq_len(burn + days)) {
    survivors <- rbinom(n_paths, occupied, survive_probability)
    arrivals <- rpois(n_paths, arrival_rate)
    demand <- survivors + arrivals
    occupied <- pmin(capacity, demand)
    if (day > burn) {
      lost_total <- lost_total + pmax(0, demand - capacity)
      arrival_total <- arrival_total + arrivals
      within_total <- within_total + (occupied <= limit)
      occupied_total <- occupied_total + occupied
    }
  }
  block_rate <- sum(lost_total) / sum(arrival_total)
  # Rejection is a ratio of patients, not the share of days with a full ward.
  # Whole paths are the independent units for the intervals; days within a
  # path are correlated because yesterday's patients may still be present.
  block_se <- sd(lost_total - block_rate * arrival_total) /
    (sqrt(n_paths) * mean(arrival_total))
  occupancy_path_prob <- within_total / days
  occupancy_prob <- mean(occupancy_path_prob)
  occupancy_se <- sd(occupancy_path_prob) / sqrt(n_paths)
  data.frame(
    beds = capacity, paths = n_paths, burn_days = burn, study_days = days,
    patient_rejection_probability = block_rate,
    patient_rejection_ci_lower = max(0, block_rate - qnorm(0.975) * block_se),
    patient_rejection_ci_upper = min(1, block_rate + qnorm(0.975) * block_se),
    probability_occupancy_at_most_85_percent = occupancy_prob,
    occupancy_ci_lower = max(0, occupancy_prob - qnorm(0.975) * occupancy_se),
    occupancy_ci_upper = min(1, occupancy_prob + qnorm(0.975) * occupancy_se),
    mean_daily_occupancy_rate = sum(occupied_total) / (n_paths * days * capacity)
  )
}

# The stationary calculation is a separate check on the simulated threshold.
# State old means the number of occupied beds after the previous admissions.
exact_hospital <- function(capacity) {
  transition <- matrix(0, capacity + 1L, capacity + 1L)
  poisson <- dpois(0:(capacity - 1L), arrival_rate)
  for (old in 0:capacity) {
    survivors <- 0:old
    survivor_prob <- dbinom(survivors, old, survive_probability)
    for (s in survivors) {
      if (s < capacity) {
        possible <- s:(capacity - 1L)
        transition[old + 1L, possible + 1L] <-
          transition[old + 1L, possible + 1L] +
          survivor_prob[s + 1L] * poisson[possible - s + 1L]
      }
    }
    # All arrival counts that fill or overflow the ward land in the full state.
    transition[old + 1L, capacity + 1L] <-
      max(0, 1 - sum(transition[old + 1L, 1:capacity]))
  }
  coefficients <- t(transition) - diag(capacity + 1L)
  coefficients[capacity + 1L, ] <- 1
  right <- c(rep(0, capacity), 1)
  # Replace one balance equation by the condition that probabilities sum to one.
  stationary <- as.vector(solve(coefficients, right))
  if (any(stationary < -1e-10)) stop("Invalid stationary distribution")
  probability <- sum(stationary[1:(floor(0.85 * capacity + 1e-10) + 1L)])
  expected_lost <- 0
  for (old in 0:capacity) {
    survivors <- 0:old
    free <- capacity - survivors
    # For Poisson arrivals A, E[(A - free)+] comes from two upper tails.
    expected_lost_conditional <- arrival_rate *
      ppois(free - 1L, arrival_rate, lower.tail = FALSE) -
      free * ppois(free, arrival_rate, lower.tail = FALSE)
    expected_lost <- expected_lost + stationary[old + 1L] *
      sum(dbinom(survivors, old, survive_probability) *
            expected_lost_conditional)
  }
  c(occupancy_probability = probability,
    patient_rejection_probability = expected_lost / arrival_rate)
}

capacity_grid <- data.frame(beds = 1:100)
exact_results <- t(vapply(capacity_grid$beds, exact_hospital, numeric(2)))
capacity_grid$exact_probability <- exact_results[, 1]
capacity_grid$exact_rejection_probability <- exact_results[, 2]
# Search by simulation on 40:70. The exact grid also checks smaller capacities,
# so an unexpectedly good result below 40 cannot be missed.
if (any(capacity_grid$exact_probability[capacity_grid$beds < 40] >= 0.95))
  stop("A capacity below the simulated grid meets the condition")
capacity_sim <- do.call(rbind, lapply(40:70, simulate_hospital))
capacity_sim$exact_probability <- capacity_grid$exact_probability[
  match(capacity_sim$beds, capacity_grid$beds)]
capacity_sim$exact_rejection_probability <- capacity_grid$exact_rejection_probability[
  match(capacity_sim$beds, capacity_grid$beds)]
eligible <- capacity_sim$beds[
  capacity_sim$probability_occupancy_at_most_85_percent >= 0.95]
if (length(eligible) == 0L) stop("No simulated capacity meets the condition")
minimum_beds <- min(eligible)
exact_minimum <- min(capacity_grid$beds[capacity_grid$exact_probability >= 0.95])
if (minimum_beds != exact_minimum) stop("Simulation and exact capacity differ")
# The decision should not depend on Monte Carlo error at 63 or 64 beds.
if (capacity_sim$occupancy_ci_lower[capacity_sim$beds == minimum_beds] < 0.95 ||
    capacity_sim$occupancy_ci_upper[capacity_sim$beds == minimum_beds - 1L] >= 0.95)
  stop("Simulation uncertainty straddles the capacity decision")
maximum_grid_difference <- max(abs(
  capacity_sim$probability_occupancy_at_most_85_percent -
    capacity_sim$exact_probability))
write.csv(data.frame(metric = c("minimum_beds_simulation", "minimum_beds_exact",
                               "maximum_absolute_grid_difference",
                               "offered_load_beds",
                               "accepted_load_beds_at_40"),
                     value = c(minimum_beds, exact_minimum,
                               maximum_grid_difference,
                               arrival_rate * 2.9,
                               arrival_rate *
                                 (1 - capacity_sim$patient_rejection_probability[
                                   capacity_sim$beds == 40]) * 2.9)),
          file.path(script_dir, "hospital_summary.csv"), row.names = FALSE)
write.csv(capacity_grid, file.path(script_dir, "hospital_capacity_exact_check.csv"),
          row.names = FALSE)
write.csv(capacity_sim, file.path(script_dir, "hospital_results.csv"), row.names = FALSE)

# This displayed year uses a separate seed and warm-up. It does not feed the
# capacity estimates above.
set.seed(1540)
occupied <- 0L
annual <- data.frame(day = seq_len(365), occupied = integer(365), arrivals = integer(365),
                     rejected = integer(365))
for (day in seq_len(730)) {
  survivors <- rbinom(1L, occupied, survive_probability)
  arrivals <- rpois(1L, arrival_rate)
  occupied <- min(40L, survivors + arrivals)
  if (day > 365) {
    index <- day - 365L
    annual$occupied[index] <- occupied
    annual$arrivals[index] <- arrivals
    annual$rejected[index] <- max(0L, survivors + arrivals - 40L)
  }
}
annual$occupancy_rate <- annual$occupied / 40
write.csv(annual, file.path(script_dir, "hospital_one_year_path.csv"), row.names = FALSE)
png(file.path(fig_dir, "hospital_one_year.png"), width = 1500, height = 900, res = 150)
plot(annual$day, annual$occupancy_rate, type = "l", col = "#245A84", lwd = 1.8,
     xlab = "Day after warm-up", ylab = "Occupancy rate (occupied beds / 40)",
     ylim = c(0.55, 1), main = "One simulated year of occupancy with 40 beds")
abline(h = 0.85, col = "#B33C3C", lty = 2, lwd = 2)
abline(h = mean(annual$occupancy_rate), col = "#2D694B", lty = 3, lwd = 2)
legend("bottomright", legend = c("85% threshold", "Year mean"),
       col = c("#B33C3C", "#2D694B"), lty = c(2, 3),
       lwd = 2, bty = "n")
dev.off()
png(file.path(fig_dir, "hospital_capacity.png"), width = 1500, height = 900, res = 150)
shown_grid <- capacity_grid[capacity_grid$beds >= 40, ]
plot(shown_grid$beds, shown_grid$exact_probability, type = "l", lwd = 1.5,
     col = "#245A84", xlab = "Number of beds",
     ylab = "Probability daily occupancy is at most 85%", ylim = c(0, 1),
     main = "Bed capacity and the occupancy requirement")
segments(capacity_sim$beds, capacity_sim$occupancy_ci_lower,
         capacity_sim$beds, capacity_sim$occupancy_ci_upper,
         col = adjustcolor("#2D694B", alpha.f = 0.65))
points(capacity_sim$beds,
       capacity_sim$probability_occupancy_at_most_85_percent,
       pch = 16, cex = 0.55, col = "#2D694B")
abline(h = 0.95, col = "#B33C3C", lty = 2, lwd = 2)
abline(v = minimum_beds, col = "#2D694B", lty = 3, lwd = 2)
legend("bottomright", c("Exact stationary", "Simulation, 95% CI",
                        "95% target"), col = c("#245A84", "#2D694B", "#B33C3C"),
       lty = c(1, NA, 2), pch = c(NA, 16, NA), bty = "n")
dev.off()

cat("Analysis completed. Minimum beds:", minimum_beds, "\n")
