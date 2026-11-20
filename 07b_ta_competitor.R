# Usage: Rscript R/07b_ta_competitor.R <dataset> <S1|S2>
# Competitor for the trait-assisted scheme (ULTE-TA; same line folds as ULTE): "UNI+secondaries".
# The target of the test lines is predicted by a regression of the target on the environment, the univariate
# genomic prediction of the target, and the observed secondary traits of the same plot.
# To avoid optimism, the genomic prediction used as covariate for the TRAINING records is itself an
# out-of-fold prediction: the training lines are split into two inner groups and each group is predicted by a
# univariate model fitted without its target phenotypes (nested cross-validation). For the test records the
# covariate is the prediction of the univariate model fitted on all training lines. In this way the covariate
# has the same nature (prediction of an unphenotyped line) in the training and test records.
# (R/07_trait_assisted.R also writes a "UNI+secondaries" row whose training covariate is the in-sample fitted
#  value; that version is optimistic in training and is NOT used in the tables.)
source("R/00_settings.R"); source("R/BMTMEMix.R")
suppressPackageStartupMessages({ library(MASS); library(nnet) })
args <- commandArgs(TRUE); d <- args[1]; scen <- args[2]
D <- readRDS(sprintf("data/prepared/%s.rds", d)); types <- scenario_types(D, scen)
eg <- eigen(D$G, symmetric = TRUE); envs <- D$envOrder
folds <- make_folds("ULTE", D$line, D$env, NFOLD, FOLD_SEED)
uni_fit <- function(Ytr, t, seed) {
  set.seed(seed)
  fu <- BMTMEMix(Ytr[, t, drop = FALSE], types[t], D$line, D$env, D$G, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                 thin = CV_MCMC$thin, eigG = eg, envLevels = envs)
  if (types[t] %in% c("nominal", "ordinal")) { up <- fu$probs_sys[[1]][fu$cell, , drop = FALSE]; as.numeric(up %*% seq_len(ncol(up))) }
  else fu$pred_sys[fu$cell, 1]
}
res <- NULL
for (f in seq_along(folds)) {
  te <- folds[[f]]$test; train <- setdiff(seq_len(nrow(D$raw)), te)
  Yd <- make_responses(D, train)
  trl <- unique(D$line[train]); set.seed(FOLD_SEED + f); ig <- sample(rep(1:2, length.out = length(trl))); names(ig) <- trl
  for (t in seq_along(types)) {
    Ytr <- Yd; Ytr[te, t] <- NA
    uni <- uni_fit(Ytr, t, FIT_SEED + 1000 * f + t)                      # test covariate
    for (k in 1:2) {                                                      # out-of-fold covariate for training records
      inner <- train[ig[D$line[train]] == k]
      Yin <- Ytr; Yin[inner, t] <- NA
      uni[inner] <- uni_fit(Yin, t, FIT_SEED + 1000 * f + 10 * k + t)[inner]
    }
    X <- data.frame(env = factor(D$env), uni = uni)
    for (k in setdiff(seq_along(types), t)) X[[names(Yd)[k]]] <- if (types[k] %in% c("nominal", "ordinal")) factor(Yd[[k]]) else Yd[[k]]
    X$y <- Yd[[t]]; Xtr <- X[train, ]; Xte <- X[te, ]
    if (types[t] %in% c("continuous", "count")) { m <- lm(y ~ ., data = Xtr); p <- predict(m, Xte); pb <- NULL }
    else if (types[t] == "binary") { m <- glm(y ~ ., data = Xtr, family = binomial("probit")); p <- predict(m, Xte, type = "response"); pb <- NULL }
    else if (types[t] == "ordinal") {
      Xtr$y <- factor(Xtr$y, ordered = TRUE)
      m <- tryCatch(polr(y ~ ., data = Xtr, method = "probit"), error = function(e) NULL)
      if (is.null(m)) m <- tryCatch(polr(y ~ ., data = Xtr, method = "logistic"), error = function(e) NULL)
      if (is.null(m)) { cat("polr failed; multinomial fallback, fold", f, "\n"); Xtr$y <- factor(Xtr$y, ordered = FALSE); m <- multinom(y ~ ., data = Xtr, trace = FALSE, maxit = 500) }
      pb <- predict(m, Xte, type = "probs"); colnames(pb) <- levels(Xtr$y); p <- NA
    } else { Xtr$y <- factor(Xtr$y); m <- multinom(y ~ ., data = Xtr, trace = FALSE, maxit = 500); pb <- predict(m, Xte, type = "probs"); colnames(pb) <- levels(Xtr$y); p <- NA }
    for (e in unique(D$env[te])) {
      w <- which(D$env[te] == e); idx <- te[w]
      mt <- metric_trait(types[t], Yd[idx, t], if (length(p) > 1) p[w] else rep(NA, length(w)), if (is.null(pb)) NULL else pb[w, , drop = FALSE])
      res <- rbind(res, data.frame(dataset = d, scenario = scen, scheme = "ULTE-TA", fold = f, env = e, trait = names(Yd)[t], type = types[t],
                                   model = "UNI+secondaries", M1 = mt[["M1"]], M2 = mt[["M2"]], M3 = mt[["M3"]], MAE = mt[["MAE"]], n = length(idx)))
    }
  }
  cat("fold", f, "done\n")
}
dir.create("outputs/cv_TA", showWarnings = FALSE, recursive = TRUE)
write.csv(res, sprintf("outputs/cv_TA/%s_%s_ULTE-TA-competitor_metrics.csv", d, scen), row.names = FALSE)
cat("done", d, scen, "\n")
