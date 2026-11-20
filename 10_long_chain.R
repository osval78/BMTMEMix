# Usage: Rscript R/10_long_chain.R <dataset>
# Stability of the cross-validation metrics with respect to the chain length: TLTE folds 1-2, scenario 2,
# BMTMEMix with the CV setting (2,500 iterations, 450 saved draws) and with chains four times longer
# (10,000 iterations, burn-in 2,800, thin 4; 1,800 saved draws), same seed.
source("R/00_settings.R"); source("R/BMTMEMix.R")
d <- commandArgs(TRUE)[1]
D <- readRDS(sprintf("data/prepared/%s.rds", d)); eg <- eigen(D$G, symmetric = TRUE); types <- D$types_S2
folds <- make_folds("TLTE", D$line, D$env, NFOLD, FOLD_SEED)[1:2]
res <- NULL
for (f in seq_along(folds)) {
  te <- folds[[f]]$test; train <- setdiff(seq_len(nrow(D$raw)), te)
  Yd <- make_responses(D, train); Ytr <- Yd; Ytr[te, ] <- NA
  for (L in list(c(2500, 700, 4), c(10000, 2800, 4))) {
    set.seed(FIT_SEED + 1000 * f); t0 <- Sys.time()
    fit <- BMTMEMix(Ytr, types, D$line, D$env, D$G, nIter = L[1], burnIn = L[2], thin = L[3], eigG = eg, envLevels = D$envOrder)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    for (t in seq_along(types)) {
      pr <- fit$pred[fit$cell, t]; pb <- if (types[t] == "ordinal") fit$probs[[t]][fit$cell, , drop = FALSE] else NULL
      for (e in unique(D$env[te])) { idx <- te[D$env[te] == e]
        mt <- metric_trait(types[t], Yd[idx, t], pr[idx], if (is.null(pb)) NULL else pb[idx, , drop = FALSE])
        res <- rbind(res, data.frame(dataset = d, fold = f, nIter = L[1], trait = names(Yd)[t], type = types[t], env = e,
                                     M1 = mt[["M1"]], M2 = mt[["M2"]], M3 = mt[["M3"]], secs = secs)) }
    }
  }
}
dir.create("outputs/long_chain", showWarnings = FALSE, recursive = TRUE)
write.csv(res, sprintf("outputs/long_chain/%s_long_chain.csv", d), row.names = FALSE)
