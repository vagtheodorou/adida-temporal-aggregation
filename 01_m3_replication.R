# ============================================================================
# 01_m3_replication.R
# Replicates Tables 2, 3, 4 of the paper (D / TD / NTD, sMAPE by bucket size),
# for both statistical and ML methods, on the M3 monthly subset.
#
# Run:  source("00_adida_core.R"); source("01_m3_replication.R")
# Output: tables printed + saved to results_m3.csv
# ============================================================================

source("00_adida_core.R")

# --- M3 lives in the Mcomp package; install from CRAN, GitHub as fallback ----
if (!requireNamespace("Mcomp", quietly = TRUE)) {
  try(install.packages("Mcomp"), silent = TRUE)
  if (!requireNamespace("Mcomp", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("robjhyndman/Mcomp")
  }
}
suppressPackageStartupMessages({
  library(Mcomp)
  library(dplyr)
  library(tidyr)
  library(gbm)
  library(xgboost)
  library(lightgbm)
  library(ranger)      # fast Random Forest (global model)
})

BUCKETS <- c(1, 2, 3, 4, 6, 12)
H       <- 12L                       # 12-month out-of-sample (paper)
USE_ID_AS_FACTOR <- FALSE            # FALSE = reproduce the paper's ML numbers
                                     # (original code used the integer id as a
                                     #  numeric feature). TRUE = the corrected,
                                     #  per-series fixed-effect version.
DROP_SEASONAL <- FALSE               # FALSE: keep seasonal series so the SD/NSD
                                     #   split is meaningful (D = all 1428 monthly).
                                     # TRUE : reproduce the paper exactly
                                     #   (D = 650 non-seasonal; then SD is empty).
M_SCALE <- 1L                        # scaling lag for scaled ME / scaled MSE:
                                     #   1  = random-walk (naive) scaling [default]
                                     #   12 = seasonal-naive scaling (monthly)

# Per-(series, method, bucket) scoring: sMAPE (accuracy), scaled ME (bias),
# scaled MSE (variation). y = training history (for the scaling denominators).
m3_score <- function(family, method, id, bucket, y, yout, f) {
  data.frame(family = family, method = method, id = id, bucket = bucket,
             smape = smape(yout, f),
             sme   = scaled_me(y, yout, f, M_SCALE),
             smse  = scaled_mse(y, yout, f, M_SCALE))
}

# ---------------------------------------------------------------------------
# Data.  subset(M3, 12) -> the 1428 MONTHLY series (Mcomp::subset.Mcomp: the
# period argument 12 == monthly, frequency 12). The ORIGINAL script then ran
# `if (SeasonalityTest(x)) next`, i.e. it DROPPED seasonal series, leaving
# D = 650 non-seasonal.  We expose that choice as DROP_SEASONAL so we can also
# keep the seasonal series and report a seasonal / non-seasonal split
# x = $x (history), xx[1:12] = test.
#
# CAVEAT: the base methods forecast the (frequency-stripped) aggregated series,
# so seasonality is NOT explicitly modelled -- exactly as in the original code.
# The SD results therefore show how temporal aggregation ALONE copes with
# seasonal series; interpret/report them in that light.
# ---------------------------------------------------------------------------
data("M3", package = "Mcomp")
M3M <- subset(M3, 12)
stopifnot(length(M3M) == 1428)       # sanity: all monthly M3 series
seasonal_all <- vapply(M3M, function(s) SeasonalityTest(as.numeric(s$x), 12), logical(1))
D   <- if (DROP_SEASONAL) M3M[!seasonal_all] else M3M
ids <- seq_along(D)

# per-series structural flags, computed on the training history $x
seasonal_flag <- vapply(D, function(s) SeasonalityTest(as.numeric(s$x), 12), logical(1))
trend_flag    <- vapply(D, function(s) is_trended(as.numeric(s$x), 0.9), logical(1))
# smoothness split (squared coefficient of variation). CV^2 = (sd/mean)^2 on the
# training history; smooth = low relative variation. 0.49 is the standard
# Syntetos-Boylan cutoff (CV^2 < 0.49 => smooth, else non-smooth/erratic).
CV2_SMOOTH_CUTOFF <- 0.49
smooth_flag <- vapply(D, function(s) {
  x <- as.numeric(s$x); m <- mean(x)
  cv2 <- if (m == 0) Inf else (sd(x) / m)^2
  is.finite(cv2) && cv2 < CV2_SMOOTH_CUTOFF
}, logical(1))
cat(sprintf("Monthly M3: %d | working set D: %d  (DROP_SEASONAL = %s)\n",
            length(M3M), length(D), DROP_SEASONAL))
cat(sprintf("  trend:  TD = %d | NTD = %d   (paper non-seasonal D: 224 / 426)\n",
            sum(trend_flag), sum(!trend_flag)))
cat(sprintf("  season: SD = %d | NSD = %d\n", sum(seasonal_flag), sum(!seasonal_flag)))
cat(sprintf("  smooth: %d | non-smooth = %d   (CV^2 < %.2f)\n",
            sum(smooth_flag), sum(!smooth_flag), CV2_SMOOTH_CUTOFF))

# ===========================================================================
# PART A — Statistical methods
# ===========================================================================
stat_rows <- list()
for (i in ids) {
  y    <- as.numeric(D[[i]]$x)
  yout <- head(as.numeric(D[[i]]$xx), H)
  for (mname in names(STAT_METHODS)) {
    fmethod <- STAT_METHODS[[mname]]
    Fmat <- matrix(NA_real_, nrow = length(BUCKETS), ncol = H,
                   dimnames = list(as.character(BUCKETS), NULL))
    for (k in BUCKETS)
      Fmat[as.character(k), ] <- safe_adida(y, H, k, fmethod)
    comb <- combine_forecasts(Fmat, exclude_k1 = TRUE)     # paper Eq. (2)
    for (k in BUCKETS)
      stat_rows[[length(stat_rows) + 1]] <-
        m3_score("stat", mname, i, as.character(k), y, yout, Fmat[as.character(k), ])
    stat_rows[[length(stat_rows) + 1]] <-
      m3_score("stat", mname, i, "Comb", y, yout, comb)
  }
}
stat_df <- bind_rows(stat_rows)

# ===========================================================================
# PART B — Global ML methods (lm, gbm, xgb, lgbm)
# Direct forecasting in aggregated space; features = lag, rolling mean, rolling sd.
# Disaggregation is equal (÷k, replicated). Combination over k = 2..12.
# ===========================================================================

# Build the global train/test design for one bucket size k.
# The last hk aggregated points of each series are the test block.
build_ml_design <- function(D, ids, k, H = 12L) {
  hk <- max(1L, round(H / k))
  all_rows <- list()
  for (i in ids) {
    y_full <- as.numeric(c(D[[i]]$x, head(as.numeric(D[[i]]$xx), H)))  # hist + test
    keep   <- (length(y_full) %/% k) * k
    AS     <- as.numeric(rollapply(tail(y_full, keep), width = k, FUN = sum, by = k))
    n      <- length(AS)
    if (n < (hk + 4)) next
    df <- data.frame(id = i, t = seq_len(n), x = AS,
                     lag = NA_real_, rmean = NA_real_, rsd = NA_real_)
    start <- max(round(20 / k) + 1, hk + 3)
    for (r in start:n) {
      df$lag[r]   <- AS[r - hk]
      df$rmean[r] <- mean(AS[(r - hk - 1):(r - hk)])
      df$rsd[r]   <- sd(AS[(r - hk - 1):(r - hk)])
    }
    df$is_test <- df$t > (n - hk)          # last hk points = test block
    all_rows[[length(all_rows) + 1]] <- df
  }
  list(data = bind_rows(all_rows), hk = hk)
}

fit_predict_ml <- function(model, train, test) {
  feats <- c("lag", "rmean", "rsd")
  if (USE_ID_AS_FACTOR) {
    train$id <- factor(train$id); test$id <- factor(test$id, levels = levels(train$id))
  }
  if (model == "lm") {
    fit <- lm(x ~ id + lag + rmean + rsd, data = na.omit(train)); predict(fit, test)
  } else if (model == "gbm") {
    fit <- gbm(x ~ id + lag + rmean + rsd, data = na.omit(train), distribution = "gaussian",
               n.trees = 100, interaction.depth = 3, verbose = FALSE)
    predict(fit, test, n.trees = 100)
  } else if (model == "xgb") {
    X  <- as.matrix(na.omit(train)[, c("id", "lag", "rmean", "rsd")])
    yv <- na.omit(train)$x
    fit <- xgb.train(params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 6,
                                   subsample = 0.8, colsample_bytree = 0.8, nthread = 2),
                     data = xgb.DMatrix(X, label = yv), nrounds = 100, verbose = 0)
    predict(fit, xgb.DMatrix(as.matrix(test[, c("id", "lag", "rmean", "rsd")])))
  } else if (model == "lgbm") {
    X  <- as.matrix(na.omit(train)[, c("id", "lag", "rmean", "rsd")])
    yv <- na.omit(train)$x
    fit <- lgb.train(params = list(objective = "regression", metric = "rmse",
                                   num_leaves = 31, learning_rate = 0.1,
                                   feature_fraction = 0.8, bagging_fraction = 0.8,
                                   bagging_freq = 5, verbose = -1),
                     data = lgb.Dataset(X, label = yv), nrounds = 100)
    predict(fit, as.matrix(test[, c("id", "lag", "rmean", "rsd")]))
  } else {  # rf -- Random Forest (global), via ranger; formula interface so
            # an id factor (when USE_ID_AS_FACTOR=TRUE) is handled natively
    fit <- ranger(x ~ id + lag + rmean + rsd, data = na.omit(train),
                  num.trees = 100, num.threads = 2, seed = 42)
    predict(fit, test)$predictions
  }
}

ml_rows <- list()
for (model in c("lm", "gbm", "xgb", "lgbm", "rf")) {
  # store per-id forecast vectors for each k, to build the combination later
  per_id_fc <- list()  # per_id_fc[[as.character(k)]][[as.character(id)]] = length-12 vector
  for (k in BUCKETS) {
    des <- build_ml_design(D, ids, k, H)
    dat <- des$data; hk <- des$hk
    train <- dat[!dat$is_test, ]; test <- dat[dat$is_test, ]
    test <- test[order(test$id, test$t), ]
    preds <- fit_predict_ml(model, train, test)
    test$frc <- pmax(preds, 0)
    # equal disaggregation: each aggregated forecast spread over k months -> 12 values
    fk <- test %>% group_by(id) %>%
      summarise(fc = list(head(rep(frc / k, each = k), H)), .groups = "drop")
    per_id_fc[[as.character(k)]] <- setNames(fk$fc, fk$id)
    for (j in seq_len(nrow(fk))) {
      i <- as.integer(fk$id[j]); y <- as.numeric(D[[i]]$x)
      yout <- head(as.numeric(D[[i]]$xx), H)
      ml_rows[[length(ml_rows) + 1]] <-
        m3_score("ml", model, i, as.character(k), y, yout, fk$fc[[j]])
    }
  }
  # combination across k = 2..12 (forecast-level average), per id
  comb_ks <- setdiff(BUCKETS, 1)
  common  <- Reduce(intersect, lapply(comb_ks, function(k) names(per_id_fc[[as.character(k)]])))
  for (idc in common) {
    mats <- sapply(comb_ks, function(k) per_id_fc[[as.character(k)]][[idc]])
    comb <- rowMeans(mats)
    i <- as.integer(idc); y <- as.numeric(D[[i]]$x)
    yout <- head(as.numeric(D[[i]]$xx), H)
    ml_rows[[length(ml_rows) + 1]] <-
      m3_score("ml", model, i, "Comb", y, yout, comb)
  }
}
ml_df <- bind_rows(ml_rows)

# ===========================================================================
# Assemble tables. Three metrics per cell, aggregated by median across series:
#   smape (accuracy) | sme (bias, signed) | smse (variation)
# ===========================================================================
all_df <- bind_rows(stat_df, ml_df)
all_df$trended  <- trend_flag[all_df$id]
all_df$seasonal <- seasonal_flag[all_df$id]
all_df$smooth   <- smooth_flag[all_df$id]

make_table <- function(df, metric = "smape") {
  df %>% group_by(method, bucket) %>%
    summarise(val = median(.data[[metric]], na.rm = TRUE), .groups = "drop") %>%
    mutate(bucket = factor(bucket, levels = c("1","2","3","4","6","12","Comb"))) %>%
    pivot_wider(names_from = bucket, values_from = val) %>%
    arrange(match(method, c("naive","naivedrift","arima","es","ces","theta",
                            "autotheta","lm","gbm","lgbm","xgb","rf")))
}

# Print sMAPE / scaled-ME / scaled-MSE tables for one subset.
METRIC_LABEL <- c(smape = "sMAPE (accuracy)",
                  sme   = "scaled ME (bias; +ve = under-forecast)",
                  smse  = "scaled MSE (variation)")
print_group <- function(df, label) {
  cat(sprintf("\n########## %s  (series = %d) ##########\n", label, length(unique(df$id))))
  if (!nrow(df)) { cat("  (empty)\n"); return(invisible()) }
  for (mt in c("smape", "sme", "smse")) {
    cat(sprintf("\n--- %s | median across series ---\n", METRIC_LABEL[[mt]]))
    print(as.data.frame(make_table(df, mt)), digits = 4)
  }
}

# Trend split (paper Tables 2-4) and Season split
print_group(all_df,                  "D  (working set)")
print_group(all_df[all_df$trended, ],   "TD  (trended)")
print_group(all_df[!all_df$trended, ],  "NTD (non-trended)")
print_group(all_df[all_df$seasonal, ],  "SD  (seasonal)")
print_group(all_df[!all_df$seasonal, ], "NSD (non-seasonal)")

# 2x2 framework (Trend x Season).
cat("\n========== 2x2 framework (Trend x Season) ==========\n")
for (tr in c(TRUE, FALSE)) for (se in c(TRUE, FALSE)) {
  sub <- all_df[all_df$trended == tr & all_df$seasonal == se, ]
  lbl <- sprintf("%s / %s",
                 ifelse(tr, "Trended", "Non-trended"),
                 ifelse(se, "Seasonal", "Non-seasonal"))
  print_group(sub, lbl)
}

write.csv(all_df, "results_m3_long.csv", row.names = FALSE)
cat("\nSaved per-series results (smape, sme, smse + trend/season/smooth flags) to results_m3_long.csv\n")
