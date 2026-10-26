# ============================================================================
# 02_m4_extension.R   (M4 Monthly + Quarterly extension)
# ----------------------------------------------------------------------------
# Runs the SAME 12 methods as the M3 replication (7 statistical + 5 ML) through
# the ADIDA sweep on M4 Monthly and Quarterly data, scored with the competition
# metrics sMAPE and MASE (plus OWA vs the naive baseline).
#
# M4 horizons are fixed by the competition: Monthly h = 18, Quarterly h = 8.
#
# Data: see load_M4() note below (M4comp2018 ships data via Git LFS, so we use
# the release tarball / direct LFS object instead of install_github).
#
# Run:  source("00_adida_core.R"); source("02_m4_extension.R")
# ============================================================================

source("00_adida_core.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
  library(gbm); library(xgboost); library(lightgbm); library(ranger)
})

# ---------------------------------------------------------------------------
# M4 data WITHOUT the broken install_github path (Git-LFS pointer problem).
# (1) GitHub RELEASE tarball, whose data/M4.rda is the real built artifact;
# (2) fallback: direct LFS object from media.githubusercontent.com.
# ---------------------------------------------------------------------------
load_M4 <- function() {
  options(timeout = max(1800, getOption("timeout")))
  ok <- requireNamespace("M4comp2018", quietly = TRUE)
  if (!ok) {
    tarball <- "https://github.com/carlanetto/M4comp2018/releases/download/0.2.0/M4comp2018_0.2.0.tar.gz"
    try(install.packages(tarball, repos = NULL, type = "source"), silent = TRUE)
    ok <- requireNamespace("M4comp2018", quietly = TRUE)
  }
  if (ok) {
    e <- new.env(); utils::data("M4", package = "M4comp2018", envir = e)
    if (exists("M4", envir = e)) return(get("M4", envir = e))
  }
  dest <- "M4.rda"
  if (!file.exists(dest) || file.size(dest) < 1e6) {
    for (u in c("https://media.githubusercontent.com/media/carlanetto/M4comp2018/master/data/M4.rda",
                "https://media.githubusercontent.com/media/carlanetto/M4comp2018/main/data/M4.rda")) {
      try(download.file(u, dest, mode = "wb"), silent = TRUE)
      if (file.exists(dest) && file.size(dest) > 1e6) break
    }
  }
  if (file.exists(dest) && file.size(dest) > 1e6) {
    e <- new.env(); load(dest, envir = e); return(get("M4", envir = e))
  }
  stop("M4 data could not be obtained. Manual route:\n",
       "  git lfs install && git clone https://github.com/carlanetto/M4comp2018.git && R CMD INSTALL M4comp2018")
}
M4 <- load_M4()
cat(sprintf("M4 loaded: %d series total\n", length(M4)))

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# h is fixed by the competition. buckets are EDITABLE (see header for options):
#   Monthly  h = 18 -> divisors-of-18 {1,2,3,6,9,18}  (alt seasonal {1,2,3,4,6,12})
#   Quarterly h = 8 -> divisors-of-8  {1,2,4,8}        (alt seasonal {1,2,4})
CONFIG <- list(
  Monthly   = list(h = 18L, m = 12L, buckets = c(1, 2, 3, 6, 9, 18)),
  Quarterly = list(h =  8L, m =  4L, buckets = c(1, 2, 4, 8))
)
SAMPLE_N         <- NA     # NA = full M4 (48k Monthly + 24k Quarterly); heavy
USE_ID_AS_FACTOR <- FALSE    # see M3 script; FALSE keeps id as a numeric feature
STAT_USE <- names(STAT_METHODS)              # all 7 statistical methods
ML_USE   <- c("lm", "gbm", "xgb", "lgbm", "rf")   # all 5 ML methods

# ===========================================================================
# Global-ML helpers (same design as the M3 script: direct forecasting in
# aggregated space; features = lag, rolling mean, rolling sd; equal disaggregation)
# ===========================================================================
build_ml_design <- function(slist, ids, k, H) {
  hk <- max(1L, round(H / k))
  all_rows <- list()
  for (i in ids) {
    y_full <- as.numeric(c(slist[[i]]$x, head(as.numeric(slist[[i]]$xx), H)))
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
    df$is_test <- df$t > (n - hk)
    all_rows[[length(all_rows) + 1]] <- df
  }
  list(data = bind_rows(all_rows), hk = hk)
}

fit_predict_ml <- function(model, train, test) {
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
    X  <- as.matrix(na.omit(train)[, c("id", "lag", "rmean", "rsd")]); yv <- na.omit(train)$x
    fit <- xgb.train(params = list(objective = "reg:squarederror", eta = 0.1, max_depth = 6,
                                   subsample = 0.8, colsample_bytree = 0.8, nthread = 2),
                     data = xgb.DMatrix(X, label = yv), nrounds = 100, verbose = 0)
    predict(fit, xgb.DMatrix(as.matrix(test[, c("id", "lag", "rmean", "rsd")])))
  } else if (model == "lgbm") {
    X  <- as.matrix(na.omit(train)[, c("id", "lag", "rmean", "rsd")]); yv <- na.omit(train)$x
    fit <- lgb.train(params = list(objective = "regression", metric = "rmse", num_leaves = 31,
                                   learning_rate = 0.1, feature_fraction = 0.8,
                                   bagging_fraction = 0.8, bagging_freq = 5, verbose = -1),
                     data = lgb.Dataset(X, label = yv), nrounds = 100)
    predict(fit, as.matrix(test[, c("id", "lag", "rmean", "rsd")]))
  } else {  # rf
    fit <- ranger(x ~ id + lag + rmean + rsd, data = na.omit(train),
                  num.trees = 100, num.threads = 2, seed = 42)
    predict(fit, test)$predictions
  }
}

# ===========================================================================
# Run one period (statistical + ML), returning a long df with sMAPE and MASE
# ===========================================================================
run_period <- function(period) {
  cfg <- CONFIG[[period]]; H <- cfg$h; m <- cfg$m; buckets <- cfg$buckets
  series <- Filter(function(s) s$period == period, M4)
  if (!is.na(SAMPLE_N) && length(series) > SAMPLE_N) {
    set.seed(1); series <- series[sample(length(series), SAMPLE_N)]
  }
  ids <- seq_along(series)
  cat(sprintf("[M4 %s] %d series, h = %d, buckets = {%s}\n",
              period, length(series), H, paste(buckets, collapse = ",")))

  score <- function(family, method, id, f) {
    y <- as.numeric(series[[id]]$x); yout <- head(as.numeric(series[[id]]$xx), H)
    data.frame(period = period, family = family, method = method, id = id,
               smape = smape(yout, f),
               mase  = if (length(y) > m) mase(y, yout, f, m) else NA_real_)
  }
  rows <- list()

  ## --- statistical (per series) -------------------------------------------
  for (i in ids) {
    y <- as.numeric(series[[i]]$x)
    for (mname in STAT_USE) {
      fmethod <- STAT_METHODS[[mname]]
      Fmat <- matrix(NA_real_, nrow = length(buckets), ncol = H,
                     dimnames = list(as.character(buckets), NULL))
      for (k in buckets) Fmat[as.character(k), ] <- safe_adida(y, H, k, fmethod)
      comb <- combine_forecasts(Fmat, exclude_k1 = TRUE)
      for (k in buckets) {
        rr <- score("stat", mname, i, Fmat[as.character(k), ]); rr$bucket <- as.character(k)
        rows[[length(rows) + 1]] <- rr
      }
      rr <- score("stat", mname, i, comb); rr$bucket <- "Comb"; rows[[length(rows) + 1]] <- rr
    }
    if (i %% 500 == 0) cat(sprintf("  stat ... %d/%d\n", i, length(series)))
  }

  ## --- ML (global per period) ---------------------------------------------
  for (model in ML_USE) {
    per_id_fc <- list()
    for (k in buckets) {
      des <- build_ml_design(series, ids, k, H)
      dat <- des$data
      train <- dat[!dat$is_test, ]; test <- dat[dat$is_test, ]
      test <- test[order(test$id, test$t), ]
      test$frc <- pmax(fit_predict_ml(model, train, test), 0)
      fk <- test %>% group_by(id) %>%
        summarise(fc = list(head(rep(frc / k, each = k), H)), .groups = "drop")
      per_id_fc[[as.character(k)]] <- setNames(fk$fc, fk$id)
      for (j in seq_len(nrow(fk))) {
        rr <- score("ml", model, as.integer(fk$id[j]), fk$fc[[j]])
        rr$bucket <- as.character(k); rows[[length(rows) + 1]] <- rr
      }
    }
    comb_ks <- setdiff(buckets, 1)
    common  <- Reduce(intersect, lapply(comb_ks, function(k) names(per_id_fc[[as.character(k)]])))
    for (idc in common) {
      comb <- rowMeans(sapply(comb_ks, function(k) per_id_fc[[as.character(k)]][[idc]]))
      rr <- score("ml", model, as.integer(idc), comb); rr$bucket <- "Comb"
      rows[[length(rows) + 1]] <- rr
    }
    cat(sprintf("  ml[%s] done\n", model))
  }
  bind_rows(rows)
}

res <- bind_rows(lapply(names(CONFIG), run_period))

# ===========================================================================
# Summary tables (competition convention: MEAN sMAPE and MASE across series)
# ===========================================================================
summ <- res %>% group_by(period, method, bucket) %>%
  summarise(sMAPE = mean(smape, na.rm = TRUE),
            MASE  = mean(mase,  na.rm = TRUE), .groups = "drop")

# OWA vs the naive baseline (method "naive", bucket "1") per period
baseline <- summ %>% filter(method == "naive", bucket == "1") %>%
  select(period, bsMAPE = sMAPE, bMASE = MASE)
summ <- summ %>% left_join(baseline, by = "period") %>%
  mutate(OWA = 0.5 * (sMAPE / bsMAPE + MASE / bMASE)) %>% select(-bsMAPE, -bMASE)

ALL_BUCKETS  <- sort(unique(unlist(lapply(CONFIG, `[[`, "buckets"))))
BUCKET_LEVELS <- c(as.character(ALL_BUCKETS), "Comb")
METHOD_ORDER  <- c("naive","naivedrift","arima","es","ces","theta","autotheta",
                   "lm","gbm","lgbm","xgb","rf")

show_metric <- function(metric) {
  summ %>% select(period, method, bucket, all_of(metric)) %>%
    mutate(bucket = factor(bucket, levels = BUCKET_LEVELS)) %>%
    pivot_wider(names_from = bucket, values_from = all_of(metric)) %>%
    arrange(period, match(method, METHOD_ORDER))
}

cat("\n=== M4 sMAPE by bucket ===\n"); print(as.data.frame(show_metric("sMAPE")), digits = 4)
cat("\n=== M4 MASE by bucket ===\n");  print(as.data.frame(show_metric("MASE")),  digits = 4)
cat("\n=== M4 OWA by bucket ===\n");   print(as.data.frame(show_metric("OWA")),   digits = 4)

write.csv(res, "results_m4_long.csv", row.names = FALSE)
cat("\nSaved per-series results to results_m4_long.csv\n")
