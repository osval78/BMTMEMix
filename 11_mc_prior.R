# =============================================================================
# 11_mc_prior.R — Monte Carlo precision of hidden-record predictions and matched marginal priors (v4)
# Usage: Rscript R/11_mc_prior.R <dataset> <S1|S2> <ULTE|TLUE|TA>
#  For fixed partitions (FOLD_SEED), every fit is repeated with three random-number streams (chains):
#  stream 1 uses the seeds of the main analysis (it reproduces version 3), streams 2 and 3 add 100000 and 200000.
#  Models: BMTMEMix with the default priors (MT), BMTMEMix with matched marginal variance priors (MT-pm;
#  priorMatch = TRUE: every variance of Sigma_A, Sigma_T and Psi has the same marginal prior as in the univariate
#  model), and the univariate models (UNI). TA = trait-assisted (target masked, secondaries of the test lines
#  observed) for the count and the categorical targets. The draws of the hidden records are saved to compute the
#  MCSE of their linear predictors.
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R"); source("R/v4_helpers.R")
args <- commandArgs(TRUE); d <- args[1]; scen <- args[2]; sch <- args[3]
D <- readRDS(sprintf("data/prepared/%s.rds", d)); types <- scenario_types(D, scen)
eg <- eigen(D$G, symmetric = TRUE); envs <- D$envOrder
covE <- covE_for(if (sch == "TA") "ULTE" else sch)
folds <- make_folds(if (sch == "TA") "ULTE" else sch, D$line, D$env, NFOLD, FOLD_SEED, envs)
targets <- if (sch == "TA") which(types %in% c("count", "nominal", "ordinal")) else NA
OUT <- list()
grab <- function(fit, k, t, te, Yd, model, stream, f, target = NA) {
  pb <- if (types[t] %in% c("nominal", "ordinal")) fit$probs[[k]][fit$cell, , drop = FALSE][te, , drop = FALSE] else NULL
  pos <- match(te, fit$draws$records)
  mc <- mc_eta(fit, k)[pos]
  dt <- data.table(dataset = d, scenario = scen, scheme = sch, stream = stream, fold = f, model = model, target = target,
                   trait = names(Yd)[t], type = types[t], line = D$line[te], env = D$env[te], obs = Yd[te, t],
                   pred = fit$pred[fit$cell, k][te], pred_sys = fit$pred_sys[fit$cell, k][te],
                   p1 = NA_real_, p2 = NA_real_, p3 = NA_real_, sd_eta = mc$sd_eta, mcse_eta = mc$mcse_eta, ess_eta = mc$ess_eta)
  if (!is.null(pb)) { dt$p1 <- pb[, 1]; dt$p2 <- pb[, 2]; dt$p3 <- pb[, 3] }
  dt
}
for (f in seq_along(folds)) {
  te <- folds[[f]]$test; mask <- union(te, folds[[f]]$drop); train <- setdiff(seq_len(nrow(D$raw)), mask)
  Yd <- make_responses(D, train)
  for (stream in 1:3) {
    off <- 100000 * (stream - 1)
    if (sch != "TA") {
      Ytr <- Yd; Ytr[mask, ] <- NA
      for (pm in c(FALSE, TRUE)) {
        set.seed(FIT_SEED + 1000 * f + off); t0 <- Sys.time()
        fit <- BMTMEMix(Ytr, types, D$line, D$env, D$G, covT = "UN", covE = covE, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                        thin = CV_MCMC$thin, eigG = eg, envLevels = envs, priorMatch = pm, drawRecords = te)
        for (t in seq_along(types)) OUT[[length(OUT) + 1]] <- grab(fit, t, t, te, Yd, if (pm) "MT-pm" else "MT", stream, f)
      }
      for (t in seq_along(types)) {
        set.seed(FIT_SEED + 1000 * f + t + off)
        fu <- BMTMEMix(Ytr[, t, drop = FALSE], types[t], D$line, D$env, D$G, covT = "UN", covE = covE, nIter = CV_MCMC$nIter,
                       burnIn = CV_MCMC$burnIn, thin = CV_MCMC$thin, eigG = eg, envLevels = envs, drawRecords = te)
        OUT[[length(OUT) + 1]] <- grab(fu, 1, t, te, Yd, "UNI", stream, f)
      }
    } else {
      for (t in targets) {
        Ytr <- Yd; Ytr[te, t] <- NA
        for (pm in c(FALSE, TRUE)) {
          set.seed(FIT_SEED + 1000 * f + 100 + t + off)
          fit <- BMTMEMix(Ytr, types, D$line, D$env, D$G, covT = "UN", covE = covE, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                          thin = CV_MCMC$thin, eigG = eg, envLevels = envs, priorMatch = pm, drawRecords = te)
          OUT[[length(OUT) + 1]] <- grab(fit, t, t, te, Yd, if (pm) "MT-pm" else "MT", stream, f, target = names(Yd)[t])
        }
        set.seed(FIT_SEED + 1000 * f + t + off)
        fu <- BMTMEMix(Ytr[, t, drop = FALSE], types[t], D$line, D$env, D$G, covT = "UN", covE = covE, nIter = CV_MCMC$nIter,
                       burnIn = CV_MCMC$burnIn, thin = CV_MCMC$thin, eigG = eg, envLevels = envs, drawRecords = te)
        OUT[[length(OUT) + 1]] <- grab(fu, 1, t, te, Yd, "UNI", stream, f, target = names(Yd)[t])
      }
    }
  }
  cat("fold", f, "done", format(Sys.time()), "\n")
}
dir.create(file.path(V4DIR, "mc_prior"), showWarnings = FALSE, recursive = TRUE)
saveRDS(rbindlist(OUT), file.path(V4DIR, "mc_prior", sprintf("%s_%s_%s.rds", d, scen, sch)))
cat("done", d, scen, sch, "\n")
