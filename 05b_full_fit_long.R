# Usage: Rscript R/05b_full_fit_long.R <dataset> <S1|S2>
# Same as 05_full_fit.R with chains five times longer (30,000 iterations, burn-in 6,000, thin 24; 1,000 draws
# per chain), to check whether the slowly mixing liability-scale variances of the discrete traits converge.
# Full-data fits: three chains with dispersed initial values and the scale move for discrete liabilities
# (pxScale = TRUE, Step 6); convergence diagnostics
# (rank-normalized split R-hat, bulk and tail ESS, Monte Carlo standard error; Vehtari et al., 2021)
# for identified quantities, selected linear predictors, and the ordinal threshold.
source("R/00_settings.R"); source("R/BMTMEMix.R")
suppressPackageStartupMessages(library(posterior))
args <- commandArgs(TRUE); d <- args[1]; scen <- args[2]
D <- readRDS(sprintf("data/prepared/%s.rds", d)); types <- scenario_types(D, scen)
eg <- eigen(D$G, symmetric = TRUE)
set.seed(FIT_SEED); tcells <- sort(sample(which(!duplicated(paste(D$line, D$env))), 4))  # 4 records to trace
fits <- list(); secs <- numeric(0)
for (ch in seq_len(FULL_MCMC$chains)) {
  set.seed(FIT_SEED + ch); t0 <- Sys.time()
  fits[[ch]] <- BMTMEMix(D$Yd, types, D$line, D$env, D$G, nIter = 30000, burnIn = 6000,
                         thin = 24, initDisp = FULL_MCMC$initDisp, eigG = eg, envLevels = D$envOrder,
                         traceCells = NULL, pxScale = TRUE)
  secs[ch] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
}
f1 <- fits[[1]]; Q <- length(f1$latnames); I <- length(f1$envs)
nm <- f1$latnames
# ---- draws of identified quantities, per chain ----
qdraws <- function(fit) {
  tr <- fit$traces; out <- list()
  SAc <- function(r) matrix(tr$SigmaA[r, ], Q, Q); SGc <- function(r) matrix(tr$SigmaG[r, ], Q, Q)
  STc <- function(r) matrix(tr$SigmaT[r, ], Q, Q); Rc <- function(r) matrix(tr$R[r, ], Q, Q)
  n <- nrow(tr$SigmaA)
  for (s in seq_len(Q)) {
    out[[paste0("SigmaA[", nm[s], "]")]] <- tr$SigmaA[, (s - 1) * Q + s]
    out[[paste0("SigmaG[", nm[s], "]")]] <- tr$SigmaG[, (s - 1) * Q + s]
    out[[paste0("SigmaTn[", nm[s], "]")]] <- tr$SigmaT[, (s - 1) * Q + s]
  }
  for (s in 1:(Q - 1)) for (l in (s + 1):Q) {
    out[[paste0("corA[", nm[s], ",", nm[l], "]")]] <- sapply(seq_len(n), function(r) cov2cor(SAc(r))[s, l])
    out[[paste0("corG[", nm[s], ",", nm[l], "]")]] <- sapply(seq_len(n), function(r) cov2cor(SGc(r))[s, l])
  }
  for (s in which(fit$types[rep(seq_along(fit$types), fit$nd)] %in% c("continuous", "count"))) out[[paste0("R[", nm[s], "]")]] <- tr$R[, (s - 1) * Q + s]
  k <- which(colSums(abs(tr$thresholds)) > 0); if (length(k)) out[["gamma_2"]] <- tr$thresholds[, k[length(k)]]
  for (s in seq_len(Q)) out[[paste0("mu[1,", nm[s], "]")]] <- tr$mu[, (s - 1) * I + 1]
  as.data.frame(out, check.names = FALSE)
}
dr <- lapply(fits, qdraws)
pars <- names(dr[[1]])
diag <- do.call(rbind, lapply(pars, function(p) {
  x <- sapply(dr, function(z) z[[p]])        # draws x chains
  data.frame(dataset = d, scenario = scen, parameter = p, mean = mean(x), sd = sd(as.vector(x)),
             q05 = unname(quantile(x, 0.05)), q95 = unname(quantile(x, 0.95)),
             rhat = rhat(x), ess_bulk = ess_bulk(x), ess_tail = ess_tail(x), mcse_mean = mcse_mean(x))
}))
# predictive quantities: linear predictors of four records and in-sample summaries of the predictions
pc <- do.call(rbind, lapply(seq_along(fits), function(ch) data.frame(chain = ch, cell = seq_len(nrow(fits[[ch]]$pred)),
         pred_cont = fits[[ch]]$pred[, 1], pred_count = fits[[ch]]$pred[, 2])))
predagree <- data.frame(dataset = d, scenario = scen,
                        cor_chains_cont = min(cor(sapply(fits, function(f) f$pred[f$cell, 1]))),
                        cor_chains_count = min(cor(sapply(fits, function(f) f$pred[f$cell, 2]))),
                        cor_chains_gv = min(cor(sapply(fits, function(f) f$gv[f$cell, 1]))),
                        max_abs_diff_prob_cat = max(abs(fits[[1]]$probs[[3]] - fits[[2]]$probs[[3]])),
                        secs_per_chain = mean(secs), draws_per_chain = nrow(dr[[1]]),
                        px_fallbacks = sum(sapply(fits, function(f) f$traces$pxFallback)))
dir.create("outputs/fullfit_long", showWarnings = FALSE, recursive = TRUE)
# pooled posterior means over chains
avg <- function(nm2) Reduce(`+`, lapply(fits, `[[`, nm2)) / length(fits)
fit <- f1
for (k in c("SigmaA", "SigmaT", "SigmaG", "SigmaE", "Omega", "R", "mu", "Ktot")) fit[[k]] <- avg(k)
fit$thresholds <- lapply(seq_along(f1$thresholds), function(t) if (is.null(f1$thresholds[[t]])) NULL else Reduce(`+`, lapply(fits, function(z) z$thresholds[[t]])) / length(fits))
fit$g <- NULL; fit$eta <- NULL
fit$draws <- dr; fit$diag <- diag; fit$predagree <- predagree
saveRDS(fit, sprintf("outputs/fullfit_long/%s_%s_fullfit.rds", d, scen))
write.csv(diag, sprintf("outputs/fullfit_long/%s_%s_diagnostics.csv", d, scen), row.names = FALSE)
write.csv(predagree, sprintf("outputs/fullfit_long/%s_%s_predagree.csv", d, scen), row.names = FALSE)
sink(sprintf("outputs/fullfit_long/%s_%s_fullfit.txt", d, scen)); print(fit)
cat("Sigma_A:\n"); print(round(fit$SigmaA, 4)); cat("Sigma_G:\n"); print(round(fit$SigmaG, 4)); cat("R:\n"); print(round(fit$R, 4))
cat("Thresholds:\n"); print(fit$thresholds); cat("Diagnostics:\n"); print(diag, digits = 3); sink()
cat("done", d, scen, "max rhat", max(diag$rhat, na.rm = TRUE), "\n")
