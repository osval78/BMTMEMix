# =============================================================================
# 15_gauss_v4.R — Gaussian multi-trait approximation with comparable treatment of uncertainty (v4)
# Usage: Rscript R/15_gauss_v4.R <dataset> <TLTE|ULTE>
#  Same partitions, seeds and chains as the Gaussian competitor of version 3 (R/03_run_cv.R ... gauss).
#  (a) "Gaussian (draw-wise)": class probabilities computed in every saved draw from the conditional predictive
#      normal of that draw, P_d(Y = c) = Phi((c + 0.5 - m_d)/s_d) - Phi((c - 0.5 - m_d)/s_d) with open extreme
#      classes, and averaged over draws; binary P(Y = 1) = mean_d Phi((m_d - 0.5)/s_d). This integrates the
#      posterior uncertainty of intercepts, genetic effects and residual variance exactly as BMTMEMix does.
#  (b) "Gaussian (plug-in)": the version-3 approximation (normal centred at the posterior mean with the posterior
#      mean residual SD).
#  (c) "Gaussian (draw-wise, training-calibrated)": the draw-wise predictive mean is recalibrated with a probit
#      (binary) or cumulative-probit (categorical) regression estimated ONLY on out-of-fold predictions of the
#      training records (inner two-fold cross-validation inside the training set; same scheme as the outer one).
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R"); source("R/v4_helpers.R")
suppressPackageStartupMessages({ library(MASS) })
args <- commandArgs(TRUE); d <- args[1]; sch <- args[2]
D <- readRDS(sprintf("data/prepared/%s.rds", d)); types <- D$types           # Gaussian model: same in S1 and S2
eg <- eigen(D$G, symmetric = TRUE); envs <- D$envOrder
folds <- make_folds(sch, D$line, D$env, NFOLD, FOLD_SEED, envs)
gfit <- function(Ytr, seed, recs) {
  set.seed(seed)
  BMTMEMix(Ytr, rep("continuous", 4), D$line, D$env, D$G, covT = "UN", covE = "UN", nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
           thin = CV_MCMC$thin, eigG = eg, envLevels = envs, drawRecords = recs)
}
gpreds <- function(fit, recs, Yd) {        # draw-wise and plug-in predictions of the binary and categorical traits
  out <- list()
  for (t in which(types %in% c("nominal", "binary"))) {
    k <- t; pos <- match(recs, fit$draws$records)
    pr <- fit$pred[fit$cell, k][recs]
    if (types[t] == "binary") {
      nd <- fit$draws$nDraw
      M <- .bt_mean(fit, k, fit$draws$condMean[pos, k, seq_len(nd)]); S <- .bt_sd(fit, k, fit$draws$condSD[pos, k, seq_len(nd)])
      out[[as.character(t)]] <- list(dw = rowMeans(pnorm((M - 0.5) / S)), plug = pmin(pmax(pr, 0), 1), mean = pr)
    } else {
      lv <- sort(unique(Yd[[t]]))
      Pd <- gauss_probs_draws(fit, k, lv)[pos, , drop = FALSE]
      sdv <- sqrt(fit$R[k, k]) * fit$csd[k]; cuts <- c(-Inf, head(lv, -1) + 0.5, Inf)
      Pp <- sapply(seq_along(lv), function(c) pnorm((cuts[c + 1] - pr) / sdv) - pnorm((cuts[c] - pr) / sdv))
      out[[as.character(t)]] <- list(dw = Pd, plug = Pp, mean = pr)
    }
  }
  out
}
OUT <- list()
for (f in seq_along(folds)) {
  te <- folds[[f]]$test; mask <- union(te, folds[[f]]$drop); train <- setdiff(seq_len(nrow(D$raw)), mask)
  Yd <- make_responses(D, train); Ytr <- Yd; Ytr[mask, ] <- NA
  fit <- gfit(Ytr, FIT_SEED + 1000 * f, te)                           # same seed as the v3 Gaussian competitor
  go <- gpreds(fit, te, Yd)
  # inner two-fold CV in the training set (records for TLTE, lines for ULTE) for the training-only recalibration
  set.seed(FOLD_SEED + 900 + f)
  if (sch == "TLTE") ig <- sample(rep(1:2, length.out = length(train))) else {
    tl <- unique(D$line[train]); g <- sample(rep(1:2, length.out = length(tl))); names(g) <- tl; ig <- g[D$line[train]] }
  oof <- lapply(go, function(z) NULL); oofm <- lapply(go, function(z) rep(NA_real_, nrow(D$raw)))
  for (k in 1:2) {
    inner <- train[ig == k]; Yin <- Ytr; Yin[inner, ] <- NA
    fi <- gfit(Yin, FIT_SEED + 1000 * f + 500 + k, inner)
    gi <- gpreds(fi, inner, Yd)
    for (nm in names(gi)) oofm[[nm]][inner] <- if (is.matrix(gi[[nm]]$dw)) as.numeric(gi[[nm]]$dw %*% seq_len(ncol(gi[[nm]]$dw))) else gi[[nm]]$mean
  }
  for (nm in names(go)) {
    t <- as.integer(nm); y <- Yd[[t]]
    xte <- if (is.matrix(go[[nm]]$dw)) as.numeric(go[[nm]]$dw %*% seq_len(ncol(go[[nm]]$dw))) else go[[nm]]$mean
    dtr <- data.frame(y = y[train], x = oofm[[nm]][train]); dte <- data.frame(x = xte)
    if (types[t] == "binary") {
      m <- glm(y ~ x, data = dtr, family = binomial("probit")); pc <- predict(m, dte, type = "response"); Pc <- NULL
    } else {
      dtr$y <- factor(dtr$y, ordered = TRUE)
      m <- tryCatch(polr(y ~ x, data = dtr, method = "probit"), error = function(e) polr(y ~ x, data = dtr, method = "logistic"))
      Pc <- predict(m, dte, type = "probs"); pc <- NA
    }
    base <- data.table(dataset = d, scheme = sch, fold = f, trait = names(Yd)[t], type = types[t], line = D$line[te], env = D$env[te], obs = y[te])
    mk <- function(model, p, P) { z <- copy(base)[, `:=`(model = model, pred = p, p1 = NA_real_, p2 = NA_real_, p3 = NA_real_)]
      if (!is.null(P)) z[, `:=`(pred = NA_real_, p1 = P[, 1], p2 = P[, 2], p3 = P[, 3])]; z }
    if (types[t] == "binary") {
      OUT[[length(OUT) + 1]] <- rbind(mk("Gaussian (plug-in)", go[[nm]]$plug, NULL), mk("Gaussian (draw-wise)", go[[nm]]$dw, NULL),
                                      mk("Gaussian (draw-wise, training-calibrated)", pc, NULL))
    } else {
      OUT[[length(OUT) + 1]] <- rbind(mk("Gaussian (plug-in)", NA, go[[nm]]$plug), mk("Gaussian (draw-wise)", NA, go[[nm]]$dw),
                                      mk("Gaussian (draw-wise, training-calibrated)", NA, Pc))
    }
  }
  cat("fold", f, "done", format(Sys.time()), "\n")
}
dir.create(file.path(V4DIR, "gauss"), showWarnings = FALSE, recursive = TRUE)
saveRDS(rbindlist(OUT), file.path(V4DIR, "gauss", sprintf("%s_%s.rds", d, sch)))
cat("done", d, sch, "\n")
