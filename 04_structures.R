# Usage: Rscript R/04_structures.R <dataset>
# Twelve covariance structures (scenario 1) evaluated on two TLTE folds with the same seed for every
# structure (common random numbers) and the same chain length as cross-validation.
#  traits      : structure of Sigma_A and Sigma_T in the main + contrast G x E model (Omega unstructured)
#  environments: structure of Sigma_E in the separable model g ~ N(0, Sigma_T x Sigma_E x G), whose
#                environmental matrix is defined on the environments themselves (trait matrix unstructured)
source("R/00_settings.R"); source("R/BMTMEMix.R")
d <- commandArgs(TRUE)[1]
D <- readRDS(sprintf("data/prepared/%s.rds", d)); eg <- eigen(D$G, symmetric = TRUE)
folds <- make_folds("TLTE", D$line, D$env, NFOLD, FOLD_SEED)[1:2]
grid <- rbind(data.frame(target = "traits", gxe = "main", covT = names(COV_CODES), covE = "UN"),
              data.frame(target = "environments", gxe = "separable", covT = "UN", covE = names(COV_CODES)))
res <- NULL
for (f in seq_along(folds)) {
  te <- folds[[f]]$test; train <- setdiff(seq_len(nrow(D$raw)), te)
  Yd <- make_responses(D, train); Ytr <- Yd; Ytr[te, ] <- NA
  for (r in seq_len(nrow(grid))) {
    set.seed(FIT_SEED + 1000 * f); t0 <- Sys.time()
    fit <- BMTMEMix(Ytr, D$types, D$line, D$env, D$G, gxe = grid$gxe[r], covT = grid$covT[r], covE = grid$covE[r],
                    nIter = STR_MCMC$nIter, burnIn = STR_MCMC$burnIn, thin = STR_MCMC$thin, eigG = eg, envLevels = D$envOrder)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    for (t in seq_along(D$types)) {
      pr <- fit$pred[fit$cell, t]; pb <- if (D$types[t] == "nominal") fit$probs[[t]][fit$cell, , drop = FALSE] else NULL
      for (e in unique(D$env[te])) {
        idx <- te[D$env[te] == e]
        mt <- metric_trait(D$types[t], Yd[idx, t], pr[idx], if (is.null(pb)) NULL else pb[idx, , drop = FALSE])
        res <- rbind(res, data.frame(dataset = d, fold = f, target = grid$target[r], gxe = grid$gxe[r], covT = grid$covT[r], covE = grid$covE[r],
                                     trait = names(Yd)[t], type = D$types[t], env = e, M1 = mt[["M1"]], M2 = mt[["M2"]], secs = secs))
      }
    }
    cat(d, f, grid$target[r], grid$covT[r], grid$covE[r], round(secs), "s\n")
  }
}
dir.create("outputs/structures", showWarnings = FALSE, recursive = TRUE)
write.csv(res, sprintf("outputs/structures/%s_structures.csv", d), row.names = FALSE)
