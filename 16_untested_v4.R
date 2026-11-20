# =============================================================================
# 16_untested_v4.R — untested environments: quantile predictive intervals and trait-independent G x E (v4)
# Usage: Rscript R/16_untested_v4.R <dataset> <TLUE|CHRONO>
#  Scenario 2, Omega = I (exchangeable), same folds and seeds as version 3.
#   MT        : BMTMEMix, unstructured Sigma_A and Sigma_T (main analysis; identical predictions to v3)
#   MT-GxEdiag: BMTMEMix with unstructured Sigma_A but diagonal Sigma_T, i.e. the G x E deviations of different
#               traits are independent, so that no common trait correlation of the G x E deviations is transferred
#               to the new environment (only the correlations of the main effects are used)
#   UNI       : univariate models
#  For continuous and count traits, 90% predictive intervals are the 5% and 95% quantiles of the posterior
#  predictive mixture (computed from the saved draws), not mean +/- 1.645 SD.
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R"); source("R/v4_helpers.R")
args <- commandArgs(TRUE); d <- args[1]; sch <- args[2]
D <- readRDS(sprintf("data/prepared/%s.rds", d)); types <- D$types_S2
eg <- eigen(D$G, symmetric = TRUE); envs <- D$envOrder
folds <- make_folds(sch, D$line, D$env, NFOLD, FOLD_SEED, envs)
OUT <- list()
grab <- function(fit, k, t, te, Yd, model, f) {
  dt <- data.table(dataset = d, scheme = sch, fold = f, model = model, trait = names(Yd)[t], type = types[t], line = D$line[te], env = D$env[te],
                   obs = Yd[te, t], pred = fit$pred[fit$cell, k][te], pred_sd = fit$pred_sd[fit$cell, k][te], q5 = NA_real_, q95 = NA_real_,
                   p1 = NA_real_, p2 = NA_real_, p3 = NA_real_)
  if (types[t] %in% c("continuous", "count")) { q <- pred_quantiles(fit, k)[match(te, fit$draws$records), , drop = FALSE]; dt$q5 <- q[, 1]; dt$q95 <- q[, 2] }
  if (types[t] == "ordinal") { pb <- fit$probs[[k]][fit$cell, , drop = FALSE][te, , drop = FALSE]; dt$p1 <- pb[, 1]; dt$p2 <- pb[, 2]; dt$p3 <- pb[, 3] }
  dt
}
for (f in seq_along(folds)) {
  te <- folds[[f]]$test; mask <- union(te, folds[[f]]$drop); train <- setdiff(seq_len(nrow(D$raw)), mask)
  Yd <- make_responses(D, train); Ytr <- Yd; Ytr[mask, ] <- NA
  for (gx in c("UN", "DIAG")) {
    set.seed(FIT_SEED + 1000 * f)
    fit <- BMTMEMix(Ytr, types, D$line, D$env, D$G, covT = "UN", covTgxe = gx, covE = "ID", nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                    thin = CV_MCMC$thin, eigG = eg, envLevels = envs, drawRecords = te)
    for (t in seq_along(types)) OUT[[length(OUT) + 1]] <- grab(fit, t, t, te, Yd, if (gx == "UN") "MT" else "MT-GxEdiag", f)
  }
  for (t in seq_along(types)) {
    set.seed(FIT_SEED + 1000 * f + t)
    fu <- BMTMEMix(Ytr[, t, drop = FALSE], types[t], D$line, D$env, D$G, covT = "UN", covE = "ID", nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                   thin = CV_MCMC$thin, eigG = eg, envLevels = envs, drawRecords = te)
    OUT[[length(OUT) + 1]] <- grab(fu, 1, t, te, Yd, "UNI", f)
  }
  cat("fold", f, "done", format(Sys.time()), "\n")
}
dir.create(file.path(V4DIR, "untested"), showWarnings = FALSE, recursive = TRUE)
saveRDS(rbindlist(OUT), file.path(V4DIR, "untested", sprintf("%s_%s.rds", d, sch)))
cat("done", d, sch, "\n")
