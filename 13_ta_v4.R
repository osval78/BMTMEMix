# =============================================================================
# 13_ta_v4.R — trait-assisted prediction of untested lines with record-level predictions (v4)
# Usage: Rscript R/13_ta_v4.R <dataset> <S1|S2> [xenv]
#  Same line folds and seeds as the ULTE-TA analysis of version 3 (R/07), so the BMTMEMix predictions are identical,
#  but every prediction is saved (fold, line, environment, model) so that paired cluster-bootstrap intervals can be
#  computed. Models:
#   MT-TA phenotypic / systematic / genetic value : BMTMEMix with only the target of the test lines masked
#   UNI                                           : univariate model of the target (target of the test lines masked)
#   UNI+sec (same plot)                           : regression of the target on the environment, the univariate genomic
#                                                   prediction and the secondary traits of the same plot (as in v3)
#   UNI+sec (across env, ridge)                   : ridge regression (glmnet, alpha = 0, lambda by line-grouped CV in the
#                                                   training set) on the environment, the genomic prediction, the same-plot
#                                                   secondaries, the means of every secondary of the line over all
#                                                   environments, and the interactions of the environment with all of them
#  In both competitors the genomic covariate of the training records is an out-of-fold prediction from an inner
#  FIVE-fold cross-validation over the training lines (each inner model uses 80% of the training lines; v3 used two
#  folds), and that of the test records the prediction of the univariate model fitted on all training lines.
#  In S2 the univariate model and the competitors are fitted only for the categorical (ordinal) target; for the
#  other targets they are identical to S1 (same data, seeds and secondaries) and are reused by R/20_summarize_v4.R.
#  xenv: the secondary traits of the test lines are also masked in the evaluated environment (BMTMEMix only).
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R"); source("R/v4_helpers.R")
suppressPackageStartupMessages({ library(MASS); library(nnet); library(glmnet) })
args <- commandArgs(TRUE); d <- args[1]; scen <- args[2]; xenv <- length(args) > 2 && args[3] == "xenv"
D <- readRDS(sprintf("data/prepared/%s.rds", d)); types <- scenario_types(D, scen)
eg <- eigen(D$G, symmetric = TRUE); envs <- D$envOrder
folds <- make_folds("ULTE", D$line, D$env, NFOLD, FOLD_SEED)
OUT <- list()
rec <- function(model, f, t, idx, Yd, pred, pb = NULL, ev = NA) {
  dt <- data.table(dataset = d, scenario = scen, scheme = if (xenv) "ULTE-TA-xenv" else "ULTE-TA", fold = f, model = model,
                   trait = names(Yd)[t], type = types[t], line = D$line[idx], env = D$env[idx], obs = Yd[idx, t], pred = pred,
                   p1 = NA_real_, p2 = NA_real_, p3 = NA_real_)
  if (!is.null(pb)) { dt$p1 <- pb[, 1]; dt$p2 <- pb[, 2]; dt$p3 <- pb[, 3] }
  dt
}
uni_fit <- function(Ytr, t, seed) {
  set.seed(seed)
  fu <- BMTMEMix(Ytr[, t, drop = FALSE], types[t], D$line, D$env, D$G, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                 thin = CV_MCMC$thin, eigG = eg, envLevels = envs)
  pb <- if (types[t] %in% c("nominal", "ordinal")) fu$probs_sys[[1]][fu$cell, , drop = FALSE] else NULL
  list(pred = if (is.null(pb)) fu$pred_sys[fu$cell, 1] else as.numeric(pb %*% seq_len(ncol(pb))), pb = pb, praw = fu$pred[fu$cell, 1])
}
fit_regression <- function(X, y, type, train, te, ridge = FALSE, foldid = NULL) {
  Xtr <- X[train, , drop = FALSE]; Xte <- X[te, , drop = FALSE]; ytr <- y[train]
  if (!ridge) {
    Xtr$y <- ytr
    if (type %in% c("continuous", "count")) { m <- lm(y ~ ., data = Xtr); return(list(p = predict(m, Xte), pb = NULL)) }
    if (type == "binary") { m <- glm(y ~ ., data = Xtr, family = binomial("probit")); return(list(p = predict(m, Xte, type = "response"), pb = NULL)) }
    if (type == "ordinal") {
      Xtr$y <- factor(Xtr$y, ordered = TRUE)
      m <- tryCatch(polr(y ~ ., data = Xtr, method = "probit"), error = function(e) NULL)
      if (is.null(m)) m <- tryCatch(polr(y ~ ., data = Xtr, method = "logistic"), error = function(e) NULL)
      if (is.null(m)) { Xtr$y <- factor(Xtr$y, ordered = FALSE); m <- multinom(y ~ ., data = Xtr, trace = FALSE, maxit = 500) }
    } else { Xtr$y <- factor(Xtr$y); m <- multinom(y ~ ., data = Xtr, trace = FALSE, maxit = 500) }
    pb <- predict(m, Xte, type = "probs"); return(list(p = NA, pb = pb))
  }
  mm <- model.matrix(~ ., data = X)[, -1, drop = FALSE]
  fam <- if (type %in% c("continuous", "count")) "gaussian" else if (type == "binary") "binomial" else "multinomial"
  cv <- cv.glmnet(mm[train, , drop = FALSE], if (fam == "multinomial") factor(ytr) else ytr, family = fam, alpha = 0, foldid = foldid,
                  type.multinomial = "grouped")
  if (fam == "multinomial") { pb <- predict(cv, mm[te, , drop = FALSE], s = "lambda.min", type = "response")[, , 1]; return(list(p = NA, pb = pb)) }
  list(p = as.numeric(predict(cv, mm[te, , drop = FALSE], s = "lambda.min", type = "response")), pb = NULL)
}
for (f in seq_along(folds)) {
  te <- folds[[f]]$test; train <- setdiff(seq_len(nrow(D$raw)), te)
  Yd <- make_responses(D, train)
  trl <- unique(D$line[train]); set.seed(FOLD_SEED + 50 + f); ig <- sample(rep(1:5, length.out = length(trl))); names(ig) <- trl
  for (t in seq_along(types)) {
    for (ev in (if (xenv) envs else "all")) {
      Ytr <- Yd; Ytr[te, t] <- NA
      if (xenv) Ytr[te[D$env[te] == ev], ] <- NA
      set.seed(FIT_SEED + 1000 * f + 100 + t + (if (xenv) 10 * match(ev, envs) else 0))
      fit <- BMTMEMix(Ytr, types, D$line, D$env, D$G, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn, thin = CV_MCMC$thin, eigG = eg, envLevels = envs)
      idx <- if (xenv) te[D$env[te] == ev] else te
      kp <- if (types[t] %in% c("nominal", "ordinal")) fit$probs[[t]][fit$cell, , drop = FALSE][idx, , drop = FALSE] else NULL
      ks <- if (types[t] %in% c("nominal", "ordinal")) fit$probs_sys[[t]][fit$cell, , drop = FALSE][idx, , drop = FALSE] else NULL
      OUT[[length(OUT) + 1]] <- rec("MT-TA phenotypic", f, t, idx, Yd, fit$pred[fit$cell, t][idx], kp)
      OUT[[length(OUT) + 1]] <- rec("MT-TA systematic", f, t, idx, Yd, fit$pred_sys[fit$cell, t][idx], ks)
      if (types[t] %in% c("continuous", "count")) OUT[[length(OUT) + 1]] <- rec("MT-TA genetic value", f, t, idx, Yd, fit$gv[fit$cell, t][idx])
    }
    if (xenv || (scen == "S2" && !types[t] %in% c("nominal", "ordinal"))) next
    # ---- univariate model and competitors ----
    u0 <- uni_fit(Ytr, t, FIT_SEED + 1000 * f + t)
    OUT[[length(OUT) + 1]] <- rec("UNI", f, t, te, Yd, if (is.null(u0$pb)) u0$praw[te] else NA, if (is.null(u0$pb)) NULL else u0$pb[te, , drop = FALSE])
    uni <- u0$pred
    for (k in 1:5) {
      inner <- train[ig[D$line[train]] == k]; Yin <- Ytr; Yin[inner, t] <- NA
      uni[inner] <- uni_fit(Yin, t, FIT_SEED + 1000 * f + 20 * k + t)$pred[inner]
    }
    sec <- setdiff(seq_along(types), t)
    X <- data.frame(env = factor(D$env), uni = uni)
    for (k in sec) X[[names(Yd)[k]]] <- if (types[k] %in% c("nominal", "ordinal")) factor(Yd[[k]]) else Yd[[k]]
    r1 <- fit_regression(X, Yd[[t]], types[t], train, te)
    OUT[[length(OUT) + 1]] <- rec("UNI+sec (same plot)", f, t, te, Yd, if (length(r1$p) > 1) r1$p else NA, r1$pb)
    # across-environment summaries of the secondaries of each line (all its observed records) + env interactions
    X2 <- data.frame(env = factor(D$env), uni = uni)
    for (k in sec) {
      v <- as.numeric(Yd[[k]]); nm <- names(Yd)[k]
      X2[[paste0(nm, "_plot")]] <- v
      X2[[paste0(nm, "_linemean")]] <- ave(v, D$line, FUN = function(z) mean(z, na.rm = TRUE))
    }
    num <- setdiff(names(X2), "env")
    for (v in num) for (e in levels(X2$env)[-1]) X2[[paste0(v, ":", e)]] <- X2[[v]] * (X2$env == e)
    fid <- ig[D$line[train]]
    r2 <- fit_regression(X2, Yd[[t]], types[t], train, te, ridge = TRUE, foldid = fid)
    OUT[[length(OUT) + 1]] <- rec("UNI+sec (across env, ridge)", f, t, te, Yd, if (length(r2$p) > 1) r2$p else NA, r2$pb)
  }
  cat("fold", f, "done", format(Sys.time()), "\n")
}
dir.create(file.path(V4DIR, "ta"), showWarnings = FALSE, recursive = TRUE)
saveRDS(rbindlist(OUT), file.path(V4DIR, "ta", sprintf("%s_%s_%s.rds", d, scen, if (xenv) "xenv" else "TA")))
cat("done", d, scen, if (xenv) "xenv" else "", "\n")
