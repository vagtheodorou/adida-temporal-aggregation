# ============================================================================
# 00_adida_core.R
# Shared engine for the ADIDA experiments (M3 replication + M4 / M5 extensions).
# Source this file first:  source("00_adida_core.R")
# ----------------------------------------------------------------------------
# Design goals: simple, reproducible, method-agnostic.
#   - ONE ADIDA function that works for any base forecaster.
#   - Base forecasters are plain functions (y, h) -> numeric vector of length h.
#   - Metrics (sMAPE, MASE, RMSSE) kept separate and explicit.
# ============================================================================

suppressPackageStartupMessages({
  library(forecast)   # rwf, auto.arima, thetaf, ses, naive, forecast()
  library(smooth)     # es, ces
  library(zoo)        # rollapply (non-overlapping aggregation)
  library(tsfeatures) # stl_features -> trend strength
})

set.seed(42)

# ---------------------------------------------------------------------------
# 1. ADIDA engine
# ---------------------------------------------------------------------------
# Aggregate -> forecast in aggregated space -> EQUAL disaggregation.
#   y           : numeric training history (original frequency)
#   h           : forecast horizon in ORIGINAL frequency
#   k           : bucket size (k = 1 means "no aggregation" = the baseline)
#   base_method : function(y_agg, h_agg) -> numeric length h_agg
# Returns a numeric vector of length h (original frequency).
adida_forecast <- function(y, h, k, base_method) {
  y <- as.numeric(y)
  if (k <= 1) {
    fk <- base_method(y, h)
    return(head(as.numeric(fk), h))
  }
  keep <- (length(y) %/% k) * k                       # drop oldest remainder
  yk   <- as.numeric(rollapply(tail(y, keep), width = k, FUN = sum, by = k))
  hk   <- max(1L, ceiling(h / k))                     # steps in aggregated space
  fk   <- as.numeric(base_method(yk, hk))
  # EQUAL disaggregation: each aggregated forecast is spread evenly over k periods
  f    <- rep(fk / k, each = k)
  head(f, h)
}

# Safe wrapper: never let one ugly series kill a 30k-series run.
safe_adida <- function(y, h, k, base_method) {
  out <- tryCatch(adida_forecast(y, h, k, base_method), error = function(e) NULL)
  if (is.null(out) || length(out) != h || any(!is.finite(out)))
    out <- rep(as.numeric(tail(y, 1)), h)             # naive fallback
  pmax(out, 0)                                        # demand is non-negative
}

# ---------------------------------------------------------------------------
# 2. Combination across bucket sizes (paper Eq. 2: average over k = 2..12)
# ---------------------------------------------------------------------------
# IMPORTANT: combine at the FORECAST level (average the forecast vectors), then
# score. This matches the paper. NOTE: the original script accidentally
# included k = 1 in the average; here we exclude it as the paper specifies.
combine_forecasts <- function(forecast_matrix, exclude_k1 = TRUE) {
  # forecast_matrix: rows = buckets (named by k as character), cols = horizon
  ks <- rownames(forecast_matrix)
  if (exclude_k1) forecast_matrix <- forecast_matrix[ks != "1", , drop = FALSE]
  colMeans(forecast_matrix)
}

# ---------------------------------------------------------------------------
# 3. Base forecasters  (each returns numeric of length h)
# ---------------------------------------------------------------------------
# All wrapped so a failure on a short aggregated series degrades to naive.
.try_fc <- function(expr, y, h) {
  out <- tryCatch(as.numeric(expr), error = function(e) NULL)
  if (is.null(out) || length(out) != h || any(!is.finite(out)))
    out <- rep(as.numeric(tail(y, 1)), h)
  out
}

bm_naive      <- function(y, h) rep(as.numeric(tail(y, 1)), h)
bm_naivedrift <- function(y, h) .try_fc(rwf(y, h = h, drift = TRUE)$mean, y, h)
bm_theta      <- function(y, h) .try_fc(thetaf(y, h = h)$mean, y, h)
bm_arima      <- function(y, h)
  .try_fc(forecast(auto.arima(y, approximation = TRUE, allowdrift = FALSE), h = h)$mean, y, h)
bm_es         <- function(y, h) .try_fc(forecast(es(y, silent = TRUE), h = h)$mean, y, h)
bm_ces        <- function(y, h) .try_fc(forecast(ces(y, silent = TRUE), h = h)$mean, y, h)
bm_autotheta  <- function(y, h) .try_fc(AutoTheta(ts(y, frequency = 1), fh = h, positive = TRUE)$mean, y, h)
bm_snaive     <- function(y, h, m = 7) .try_fc(rep(as.numeric(tail(y, m)), length.out = h), y, h)

# Registry used by the statistical experiments
STAT_METHODS <- list(
  naive      = bm_naive,
  naivedrift = bm_naivedrift,
  es         = bm_es,
  ces        = bm_ces,
  theta      = bm_theta,
  arima      = bm_arima,
  autotheta  = bm_autotheta
)

# ---------------------------------------------------------------------------
# 4. Metrics
# ---------------------------------------------------------------------------
smape <- function(a, f) {
  a <- as.numeric(a); f <- as.numeric(f)
  d <- abs(a) + abs(f)
  term <- ifelse(d == 0, 0, 2 * abs(f - a) / d)
  mean(term) * 100
}

# MASE with seasonal-naive in-sample scaling (M4 competition convention).
mase <- function(train, a, f, m) {
  scal <- mean(abs(diff(as.numeric(train), lag = m)))
  if (!is.finite(scal) || scal == 0) return(NA_real_)
  mean(abs(as.numeric(a) - as.numeric(f))) / scal
}

# Scaled Mean Error  -> BIAS. Signed: positive = under-forecasting (actual > fc).
# Denominator = in-sample naive MAE at lag m (same scaling family as MASE).
scaled_me <- function(train, a, f, m = 1) {
  scal <- mean(abs(diff(as.numeric(train), lag = m)))
  if (!is.finite(scal) || scal == 0) return(NA_real_)
  mean(as.numeric(a) - as.numeric(f)) / scal
}

# Scaled MSE -> VARIATION / squared accuracy. Denominator = in-sample naive MSE
# at lag m. Equivalent to RMSSE^2 when m matches and no zero-trimming is applied.
scaled_mse <- function(train, a, f, m = 1) {
  scal <- mean(diff(as.numeric(train), lag = m)^2)
  if (!is.finite(scal) || scal == 0) return(NA_real_)
  mean((as.numeric(a) - as.numeric(f))^2) / scal
}

# RMSSE scaling denominator (M5): mean squared 1-step diff over the training
# series. DEFAULT trim_zeros = TRUE matches the official M5 scorer, which trims
# each series to its first non-zero sale before mean(diff(.)^2). Verified by
# direct computation: evaluation+trim reproduces the published Naive L12 = 1.2534.
# trim_zeros = FALSE uses the full series. The M5 script passes SCALE_TRIM_ZEROS.
rmsse_scale <- function(train, trim_zeros = TRUE) {
  train <- as.numeric(train)
  if (trim_zeros) {
    nz <- which(train != 0)
    if (length(nz) == 0) return(NA_real_)
    train <- train[nz[1]:length(train)]
  }
  if (length(train) < 2) return(NA_real_)
  mean(diff(train)^2)
}
rmsse <- function(train, a, f, trim_zeros = TRUE) {
  s <- rmsse_scale(train, trim_zeros)
  if (!is.finite(s) || s == 0) return(NA_real_)
  sqrt(mean((as.numeric(a) - as.numeric(f))^2) / s)
}

# ---------------------------------------------------------------------------
# 5. Trend strength (for the D / TD / NTD split)
# ---------------------------------------------------------------------------
# Use NAMED extraction, not positional [3]. Positional indexing is fragile:
# tsfeatures::stl_features() column order/length changes with version and with
# whether the series is seasonal, so [3] can silently grab "linearity".
trend_strength <- function(x) {
  v <- tsfeatures::stl_features(ts(as.numeric(x), frequency = 1))
  as.numeric(v["trend"])
}
is_trended <- function(x, thresh = 0.9) {
  ts_ <- trend_strength(x)
  isTRUE(ts_ > thresh)
}

# ---------------------------------------------------------------------------
# 6. Theta family + seasonality test (ported verbatim from the original script,
#    these form the AutoTheta engine and are correct as-is).
# ---------------------------------------------------------------------------
SeasonalityTest <- function(input, ppy = 12) {
  tcrit <- 1.645
  if (length(input) < 3 * ppy) {
    test_seasonal <- FALSE
  } else {
    xacf <- acf(input, plot = FALSE)$acf[-1, 1, 1]
    clim <- tcrit / sqrt(length(input)) * sqrt(cumsum(c(1, 2 * xacf^2)))
    test_seasonal <- (abs(xacf[ppy]) > clim[ppy])
    if (is.na(test_seasonal)) test_seasonal <- FALSE
  }
  test_seasonal
}

Theta.models.fit <- function(input, fh, theta, curve, model, seasonality,
                             plot = FALSE, positive = TRUE) {
  if (theta < 1) theta <- 1
  if (fh < 1) fh <- 1
  outtest <- naive(input, h = fh)$mean
  wses <- (1 / theta); wlrl <- (1 - wses)
  ppy <- frequency(input)
  if (seasonality == "N") {
    des_input <- input; SIout <- rep(1, fh); SIin <- rep(1, length(input))
  } else if (seasonality == "A") {
    Dec <- decompose(input, type = "additive")
    des_input <- input - Dec$seasonal
    SIin <- Dec$seasonal
    SIout <- head(rep(Dec$seasonal[(length(Dec$seasonal) - ppy + 1):length(Dec$seasonal)], fh), fh)
  } else {
    Dec <- decompose(input, type = "multiplicative")
    des_input <- input / Dec$seasonal
    SIin <- Dec$seasonal
    SIout <- head(rep(Dec$seasonal[(length(Dec$seasonal) - ppy + 1):length(Dec$seasonal)], fh), fh)
  }
  observations <- length(des_input)
  xs <- c(1:observations)
  xf <- xff <- c((observations + 1):(observations + fh))
  newdf <- data.frame(xs = xff)
  if (curve == "Exp") {
    estimate <- lm(log(des_input) ~ xs)
    thetaline0In  <- exp(predict(estimate)) + input - input
    thetaline0Out <- exp(predict(estimate, newdf)) + outtest - outtest
  } else {
    estimate <- lm(des_input ~ poly(xs, 1, raw = TRUE))
    thetaline0In  <- predict(estimate) + des_input - des_input
    thetaline0Out <- predict(estimate, newdf) + outtest - outtest
  }
  if (model == "A") {
    thetalineT <- theta * des_input + (1 - theta) * thetaline0In
  } else {
    thetalineT <- (des_input^theta) * (thetaline0In^(1 - theta))
  }
  sesmodel <- ses(thetalineT, h = fh)
  thetaline2In  <- sesmodel$fitted
  thetaline2Out <- sesmodel$mean
  if (model == "A") {
    forecastsIn  <- as.numeric(thetaline2In * wses)  + as.numeric(thetaline0In * wlrl)  + des_input - des_input
    forecastsOut <- as.numeric(thetaline2Out * wses) + as.numeric(thetaline0Out * wlrl) + outtest - outtest
  } else {
    forecastsIn  <- ((as.numeric(thetaline2In)^(1 / theta))  * (as.numeric(thetaline0In)^(1 - (1 / theta))))  + des_input - des_input
    forecastsOut <- ((as.numeric(thetaline2Out)^(1 / theta)) * (as.numeric(thetaline0Out)^(1 - (1 / theta)))) + outtest - outtest
  }
  if (seasonality == "A") {
    forecastsIn  <- forecastsIn + SIin
    forecastsOut <- forecastsOut + SIout
  } else {
    forecastsIn  <- forecastsIn * SIin
    forecastsOut <- forecastsOut * SIout
  }
  if (positive) forecastsOut[forecastsOut < 0] <- 0
  list(fitted = forecastsIn, mean = forecastsOut,
       fitted0 = thetaline0In, mean0 = thetaline0Out,
       fitted2 = thetaline2In, mean2 = thetaline2Out,
       model = paste(seasonality, model, curve, round(theta, 2),
                     round(sesmodel$model$par[1], 3), round(sesmodel$model$par[2], 3)))
}

AutoTheta <- function(input, fh, positive = TRUE) {
  if (min(input) > 0) {
    molist <- c("M", "A"); trlist <- c("Lrl", "Exp")
  } else {
    molist <- c("A"); trlist <- c("Lrl")
  }
  base <- mean(input); input <- input / base
  ppy <- frequency(input); ST <- FALSE
  if (ppy > 1) ST <- SeasonalityTest(input, ppy)
  if (ST) {
    selist <- c("M", "A"); listnames <- c()
    for (i in seq_along(selist))
      for (ii in seq_along(molist))
        for (iii in seq_along(trlist))
          listnames <- c(listnames, paste(selist[i], molist[ii], trlist[iii]))
  } else {
    listnames <- c()
    for (ii in seq_along(molist))
      for (iii in seq_along(trlist))
        listnames <- c(listnames, paste("N", molist[ii], trlist[iii]))
  }
  excluded <- c("N M Lrl", "A M Lrl", "A M Exp", "M M Lrl")
  listnames <- listnames[!(listnames %in% excluded)]
  modellist <- NULL
  for (i in seq_along(listnames))
    modellist[[length(modellist) + 1]] <-
      c(substr(listnames, 1, 1)[i], substr(listnames, 3, 3)[i], substr(listnames, 5, 7)[i])
  errorsin <- c(); models <- NULL
  optfun <- function(x, input, fh, curve, model, seasonality)
    mean(abs(Theta.models.fit(input, fh, theta = x, curve, model, seasonality, plot = FALSE)$fitted - input))
  for (j in seq_along(listnames)) {
    optTheta <- suppressWarnings(optimize(optfun, c(1, 3), input = input, fh = fh,
                                          curve = modellist[[j]][3], model = modellist[[j]][2],
                                          seasonality = modellist[[j]][1])$minimum)
    fortheta <- Theta.models.fit(input, fh = fh, theta = optTheta, curve = modellist[[j]][3],
                                 model = modellist[[j]][2], seasonality = modellist[[j]][1], plot = FALSE)
    models[[length(models) + 1]] <- fortheta
    errorsin <- c(errorsin, mean(abs(input - fortheta$fitted)))
  }
  selected.model <- models[[which.min(errorsin)]]
  frc <- selected.model$mean * base
  if (positive) frc[frc < 0] <- 0
  list(mean = frc, description = selected.model$model)
}

message("00_adida_core.R loaded.")
