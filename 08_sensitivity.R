# Usage: Rscript R/08_sensitivity.R <dataset>
# Sensitivity of BMTMEMix (scenario 2, TLTE fold 1) to (i) the prior (R2 and inverse-Wishart degrees of
# freedom), and (ii) the order of traits and environments for the order-dependent structures
# (AR1, ARH1, TOEP, ANTE1, SP), comparing the original and the reversed order.
source("R/00_settings.R"); source("R/BMTMEMix.R")
d <- commandArgs(TRUE)[1]
D <- readRDS(sprintf("data/prepared/%s.rds", d)); eg <- eigen(D$G, symmetric = TRUE)
types <- D$types_S2
fo <- make_folds("TLTE", D$line, D$env, NFOLD, FOLD_SEED)[[1]]
te <- fo$test; train <- setdiff(seq_len(nrow(D$raw)), te)
Yd <- make_responses(D, train); Ytr <- Yd; Ytr[te, ] <- NA
evalfit <- function(fit, ord = seq_along(types)) {
  out <- NULL
  for (k in seq_along(ord)) {
    t <- ord[k]; pr <- fit$pred[fit$cell, k]; pb <- if (types[t] == "ordinal") fit$probs[[k]][fit$cell, , drop = FALSE] else NULL
    for (e in unique(D$env[te])) { idx <- te[D$env[te] == e]
      mt <- metric_trait(types[t], Yd[idx, t], pr[idx], if (is.null(pb)) NULL else pb[idx, , drop = FALSE])
      out <- rbind(out, data.frame(trait = names(Yd)[t], type = types[t], env = e, M1 = mt[["M1"]], M2 = mt[["M2"]], M3 = mt[["M3"]])) }
  }
  out
}
res <- NULL
prior_grid <- data.frame(R2 = c(0.5, 0.25, 0.75, 0.5, 0.5), dfIW = c(4, 4, 4, 2, 20))
for (r in seq_len(nrow(prior_grid))) {
  set.seed(FIT_SEED + 1000)
  fit <- BMTMEMix(Ytr, types, D$line, D$env, D$G, R2 = prior_grid$R2[r], dfIW = prior_grid$dfIW[r],
                  nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn, thin = CV_MCMC$thin, eigG = eg, envLevels = D$envOrder)
  ev <- evalfit(fit); ev$analysis <- "prior"; ev$setting <- sprintf("R2=%.2f, df=q+%d", prior_grid$R2[r], prior_grid$dfIW[r])
  ev$corG_12 <- cov2cor(fit$SigmaG)[1, 2]; ev$corG_23 <- cov2cor(fit$SigmaG)[2, 3]; ev$corA_34 <- cov2cor(fit$SigmaA)[3, 4]
  res <- rbind(res, ev); cat("prior", r, "\n")
}
for (st in c("AR1", "ARH1", "TOEP", "ANTE1", "SP")) for (rev in c(FALSE, TRUE)) {
  # traits (main model)
  ord <- if (rev) rev(seq_along(types)) else seq_along(types)
  set.seed(FIT_SEED + 1000)
  fit <- BMTMEMix(Ytr[, ord], types[ord], D$line, D$env, D$G, covT = st, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                  thin = CV_MCMC$thin, eigG = eg, envLevels = D$envOrder)
  ev <- evalfit(fit, ord); ev$analysis <- "order-traits"; ev$setting <- paste(st, if (rev) "reversed" else "original")
  ev$corG_12 <- NA; ev$corG_23 <- NA; ev$corA_34 <- NA; res <- rbind(res, ev)
  # environments (separable model)
  el <- if (rev) rev(D$envOrder) else D$envOrder
  set.seed(FIT_SEED + 1000)
  fit <- BMTMEMix(Ytr, types, D$line, D$env, D$G, gxe = "separable", covE = st, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                  thin = CV_MCMC$thin, eigG = eg, envLevels = el)
  ev <- evalfit(fit); ev$analysis <- "order-environments"; ev$setting <- paste(st, if (rev) "reversed" else "original")
  ev$corG_12 <- NA; ev$corG_23 <- NA; ev$corA_34 <- NA; res <- rbind(res, ev); cat("order", st, rev, "\n")
}
res$dataset <- d
dir.create("outputs/sensitivity", showWarnings = FALSE, recursive = TRUE)
write.csv(res, sprintf("outputs/sensitivity/%s_sensitivity.csv", d), row.names = FALSE)
