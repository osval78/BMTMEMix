# Usage: Rscript R/07_trait_assisted.R <dataset> <S1|S2> [xenv]
# Trait-assisted prediction of untested lines (ULTE-TA). Same line folds as ULTE.
#  default: the target trait of the test lines is masked in every environment; their other traits stay observed.
#  xenv   : additionally, the other traits of the test lines are masked in the environment being evaluated,
#           so that the evaluated phenotype shares no residual with any observed phenotype (Runcie and Cheng, 2019).
# Predictors stored for BMTMEMix: phenotypic (conditional on the other latent dimensions of the cell),
# systematic (eta, marginal residual) and genetic value. Competitor (same secondaries):
#  "UNI+secondaries": regression of the target on the environment, the univariate genomic prediction and the
#  observed secondary traits, fitted on training cells (lm / probit glm / proportional-odds or multinomial model).
#  NOTE: here the univariate covariate of the training cells is an in-sample fitted value, which makes this
#  competitor optimistic in training and weak in test; the tables use the nested-CV version of
#  R/07b_ta_competitor.R instead (these rows are dropped by R/06_summarize.R).
source("R/00_settings.R"); source("R/BMTMEMix.R")
suppressPackageStartupMessages({ library(MASS); library(nnet) })
args <- commandArgs(TRUE); d <- args[1]; scen <- args[2]; xenv <- length(args) > 2 && args[3] == "xenv"
D <- readRDS(sprintf("data/prepared/%s.rds", d)); types <- scenario_types(D, scen)
eg <- eigen(D$G, symmetric = TRUE)
folds <- make_folds("ULTE", D$line, D$env, NFOLD, FOLD_SEED)
envs <- D$envOrder
res <- NULL; preds <- NULL
metr <- function(type, y, p, pb) metric_trait(type, y, p, pb)
add <- function(model, f, e, t, idx, Yd, p, pb) {
  mt <- metr(types[t], Yd[idx, t], p, pb)
  data.frame(dataset = d, scenario = scen, scheme = if (xenv) "ULTE-TA-xenv" else "ULTE-TA", fold = f, env = e,
             trait = names(Yd)[t], type = types[t], model = model, M1 = mt[["M1"]], M2 = mt[["M2"]], M3 = mt[["M3"]],
             MAE = mt[["MAE"]], n = length(idx))
}
for (f in seq_along(folds)) {
  te <- folds[[f]]$test
  train <- setdiff(seq_len(nrow(D$raw)), te)
  Yd <- make_responses(D, train)                      # thresholds from training lines only
  for (t in seq_along(types)) {
    evals <- if (xenv) envs else "all"
    for (ev in evals) {
      Ytr <- Yd; Ytr[te, t] <- NA
      if (xenv) Ytr[te[D$env[te] == ev], ] <- NA        # no secondary of the test lines in the evaluated environment
      set.seed(FIT_SEED + 1000 * f + 100 + t + (if (xenv) 10 * match(ev, envs) else 0))
      fit <- BMTMEMix(Ytr, types, D$line, D$env, D$G, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                      thin = CV_MCMC$thin, eigG = eg, envLevels = envs)
      kp <- if (types[t] %in% c("nominal", "ordinal")) fit$probs[[t]][fit$cell, , drop = FALSE] else NULL
      ks <- if (types[t] %in% c("nominal", "ordinal")) fit$probs_sys[[t]][fit$cell, , drop = FALSE] else NULL
      for (e in (if (xenv) ev else unique(D$env[te]))) {
        idx <- te[D$env[te] == e]
        res <- rbind(res, add("BMTMEMix (phenotypic)", f, e, t, idx, Yd, fit$pred[fit$cell, t][idx], if (is.null(kp)) NULL else kp[idx, , drop = FALSE]),
                          add("BMTMEMix (systematic)", f, e, t, idx, Yd, fit$pred_sys[fit$cell, t][idx], if (is.null(ks)) NULL else ks[idx, , drop = FALSE]))
        if (types[t] %in% c("continuous", "count"))
          res <- rbind(res, add("BMTMEMix (genetic value)", f, e, t, idx, Yd, fit$gv[fit$cell, t][idx], NULL))
      }
      if (!xenv) {
        # ---- competitor: univariate genomic prediction + regression on the observed secondary traits ----
        set.seed(FIT_SEED + 1000 * f + t)
        fu <- BMTMEMix(Ytr[, t, drop = FALSE], types[t], D$line, D$env, D$G, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                       thin = CV_MCMC$thin, eigG = eg, envLevels = envs)
        uni <- fu$pred_sys[fu$cell, 1]
        if (types[t] %in% c("nominal", "ordinal")) { up <- fu$probs_sys[[1]][fu$cell, , drop = FALSE]; uni <- as.numeric(up %*% seq_len(ncol(up))) }
        X <- data.frame(env = factor(D$env), uni = uni)
        for (k in setdiff(seq_along(types), t)) X[[names(Yd)[k]]] <- if (types[k] %in% c("nominal", "ordinal")) factor(Yd[[k]]) else Yd[[k]]
        X$y <- Yd[[t]]; Xtr <- X[train, ]; Xte <- X[te, ]
        if (types[t] %in% c("continuous", "count")) { m <- lm(y ~ ., data = Xtr); p <- predict(m, Xte); pb <- NULL }
        else if (types[t] == "binary") { m <- glm(y ~ ., data = Xtr, family = binomial("probit")); p <- predict(m, Xte, type = "response"); pb <- NULL }
        else if (types[t] == "ordinal") {
          # proportional-odds probit; logistic link, then multinomial model if the probit fit fails to start/converge
          Xtr$y <- factor(Xtr$y, ordered = TRUE)
          m <- tryCatch(polr(y ~ ., data = Xtr, method = "probit"), error = function(e) NULL)
          if (is.null(m)) m <- tryCatch(polr(y ~ ., data = Xtr, method = "logistic"), error = function(e) NULL)
          if (is.null(m)) { cat("ordinal competitor: polr failed, multinomial fallback (fold", f, ")\n"); Xtr$y <- factor(Xtr$y, ordered = FALSE); m <- multinom(y ~ ., data = Xtr, trace = FALSE, maxit = 500) }
          pb <- predict(m, Xte, type = "probs"); colnames(pb) <- levels(Xtr$y); p <- NA }
        else { Xtr$y <- factor(Xtr$y); m <- multinom(y ~ ., data = Xtr, trace = FALSE, maxit = 500); pb <- predict(m, Xte, type = "probs"); colnames(pb) <- levels(Xtr$y); p <- NA }
        for (e in unique(D$env[te])) {
          w <- which(D$env[te] == e); idx <- te[w]
          res <- rbind(res, add("UNI+secondaries", f, e, t, idx, Yd, if (length(p) > 1) p[w] else rep(NA, length(w)), if (is.null(pb)) NULL else pb[w, , drop = FALSE]))
        }
      }
    }
  }
  cat("fold", f, "done\n")
}
dir.create("outputs/cv_TA", showWarnings = FALSE, recursive = TRUE)
write.csv(res, sprintf("outputs/cv_TA/%s_%s_%s_metrics.csv", d, scen, if (xenv) "ULTE-TA-xenv" else "ULTE-TA"), row.names = FALSE)
cat("done", d, scen, if (xenv) "xenv" else "", "\n")
