# =============================================================================
#  BMTMEMix — Bayesian Mixed-response Multi-Trait Multi-Environment model
#  R interface (version 3). The Gibbs sampler is implemented in C++ (src/BMTMEMix_core.cpp).
# =============================================================================
suppressPackageStartupMessages({ library(Rcpp); library(RcppArmadillo) })

.bmtme_src <- function() {
  cand <- c(file.path("src", "BMTMEMix_core.cpp"), file.path("..", "src", "BMTMEMix_core.cpp"))
  cand[file.exists(cand)][1]
}
if (!exists("bmtmemix_gibbs")) {
  .t0 <- Sys.time()
  Rcpp::sourceCpp(.bmtme_src(), cacheDir = file.path(tempdir(), "bmtme_cache"))
  BMTME_COMPILE_SECONDS <- as.numeric(difftime(Sys.time(), .t0, units = "secs"))
}

COV_CODES <- c(UN = 1, DIAG = 2, FA = 3, REC = 4, CS = 5, ID = 6, CSH = 7,
               AR1 = 8, ARH1 = 9, TOEP = 10, ANTE1 = 11, SP = 12)
TYPE_CODES <- c(continuous = 0, binary = 1, nominal = 2, count = 3, ordinal = 4)

## genomic relationship matrix (VanRaden-type, standardized markers) ---------
make_G <- function(M, ridge = 0.01) {
  M <- as.matrix(M)
  M <- M[, apply(M, 2, function(x) var(x) > 0), drop = FALSE]
  p <- colMeans(M) / 2
  M <- M[, pmin(p, 1 - p) >= 0.05, drop = FALSE]
  W <- scale(M, center = TRUE, scale = TRUE)
  G <- tcrossprod(W) / ncol(W)
  G <- G / mean(diag(G))
  G <- G + ridge * diag(nrow(G))
  attr(G, "n_used") <- ncol(W)
  G
}

#' Fit BMTMEMix
#' @param Y      data.frame of responses (rows = records, one column per trait)
#' @param types  "continuous", "binary", "nominal", "count", "ordinal" (ordinal classes ordered by their codes)
#' @param line,env  line and environment of each record
#' @param G      genomic relationship matrix (rownames = line ids)
#' @param gxe    "main" (default): main genomic effect + G x E deviations in the space of environmental
#'               contrasts (identified decomposition; covE = "UN" or "ID" = exchangeable environments);
#'               "separable": g ~ N(0, SigmaT x SigmaE x G), any of the 12 structures for SigmaE.
#' @param df0    degrees of freedom of the scaled inverse chi-square priors (BGLR default 5)
#' @param dfIW   inverse-Wishart priors use df = q + dfIW (default 4, i.e. 5 when q = 1); dfIW > 0 keeps them proper
#' @param initDisp dispersion of the initial values (0 = data-based inits; > 0 for multiple dispersed chains)
#' @param envLevels optional order of the environments (default: sorted unique values)
#' @param pxScale  TRUE adds a scale (group) move for binary and ordinal liabilities (Liu and Sabatti, 2000) that
#'               improves the mixing of their liability-scale variances and thresholds; FALSE (default) leaves the
#'               sampler and its random-number stream unchanged.
BMTMEMix <- function(Y, types, line, env, G, gxe = c("main", "separable"), covT = "UN", covE = "UN",
                     nIter = 6000, burnIn = 2000, thin = 5, R2 = 0.5, df0 = 5, dfIW = 4,
                     nfT = 1, nfE = 1, offset = NULL, DT = NULL, DE = NULL, eigG = NULL,
                     gammaMax = 20, initDisp = 0, traceCells = NULL, envLevels = NULL,
                     m0 = NULL, vm0 = NULL, priorScales = NULL, nPredDraws = 10, pxScale = FALSE, verbose = FALSE,
                     priorMatch = FALSE, fixCov = NULL, drawRecords = NULL, covTgxe = NULL) {
  #  covTgxe: structure of the G x E trait matrix Sigma_T when it differs from that of Sigma_A (e.g. "DIAG":
  #    trait-independent G x E deviations; the main-effect matrix Sigma_A keeps the structure covT).
  # v4 options:
  #  priorMatch = TRUE: the scale matrices of the unstructured inverse-Wishart priors (Sigma_A, Sigma_T, Psi) use the
  #    multiplier (dfIW + 3) of the one-dimensional case instead of (q + dfIW + q + 1), so that the marginal prior of
  #    every variance, IG(dfIW/2 + 1/2, (dfIW + 3) S0 / 2), is the same in the multi-trait and in the univariate model.
  #  fixCov = list(SigmaA, SigmaT, Omega, R): covariance matrices fixed at known values (latent internal scale;
  #    oracle fits used in the simulation of Section S8); intercepts, effects, liabilities and thresholds are sampled.
  #  drawRecords: record indices whose draw-level predictive moments (conditional and systematic means and SDs of the
  #    latent dimension of each trait, and the genetic value) are returned in fit$draws.
  gxe <- match.arg(gxe)
  Y <- as.data.frame(Y); T <- ncol(Y)
  stopifnot(length(types) == T, all(types %in% names(TYPE_CODES)))
  if (!covT %in% names(COV_CODES) || !covE %in% names(COV_CODES)) stop("unknown covariance structure")
  if (dfIW <= 0) stop("dfIW must be > 0 so that every inverse-Wishart prior is proper (df > q - 1)")
  if (df0 <= 0) stop("df0 must be > 0")
  tcode <- as.integer(TYPE_CODES[types])
  lines <- rownames(G); envs <- if (is.null(envLevels)) sort(unique(as.character(env))) else envLevels
  J <- length(lines); I <- length(envs)
  jj <- match(as.character(line), lines); ii <- match(as.character(env), envs)
  if (anyNA(jj) || anyNA(ii)) stop("some lines are not in G or some environments are not in envLevels")
  cell <- jj + J * (ii - 1)
  Ym <- matrix(NaN, J * I, T); Om <- matrix(0, J * I, T)
  for (t in seq_len(T)) { Ym[cell, t] <- as.numeric(Y[[t]]); if (!is.null(offset)) Om[cell, t] <- offset[, t] }
  Ym[is.na(Ym)] <- NaN
  Kcat <- rep(1L, T); levs <- vector("list", T); cmean <- rep(0, T); csd <- rep(1, T)
  for (t in seq_len(T)) {
    y <- Ym[, t]; o <- is.finite(y)
    if (types[t] == "continuous") { cmean[t] <- mean(y[o]); csd[t] <- sd(y[o]); Ym[o, t] <- (y[o] - cmean[t]) / csd[t] }
    else if (types[t] %in% c("nominal", "ordinal")) { levs[[t]] <- sort(unique(y[o])); Kcat[t] <- length(levs[[t]]); Ym[o, t] <- match(y[o], levs[[t]]) }
    else if (types[t] == "binary") stopifnot(all(y[o] %in% c(0, 1)))
    else if (types[t] == "count") stopifnot(all(y[o] >= 0), all(abs(y[o] - round(y[o])) < 1e-8))
  }
  envc <- rep(seq_len(I), each = J)
  pooled_var <- function(z) { o <- is.finite(z); v <- tapply(z[o], envc[o], var); mean(v[is.finite(v)]) }
  nd <- ifelse(types == "nominal", Kcat - 1, 1); Q <- sum(nd)
  vlat <- unlist(lapply(seq_len(T), function(t) {
    if (types[t] == "continuous") rep(pooled_var(Ym[, t]), nd[t])
    else if (types[t] == "count") rep(pooled_var(log(Ym[, t] + 0.5) - Om[, t]), nd[t])
    else rep(1, nd[t]) }))
  vlat[!is.finite(vlat) | vlat <= 0] <- 1
  useMain <- gxe == "main" && I > 1
  s0A <- R2 * vlat * 0.5; s0T <- R2 * vlat * (if (useMain) 0.5 else 1); s0R <- (1 - R2) * vlat
  IE <- if (useMain) I - 1 else I
  if (!is.null(priorScales)) { s0A <- priorScales$s0A; s0T <- priorScales$s0T; s0R <- priorScales$s0R }   # explicit prior scales (e.g. SBC)
  nfr <- sum(types %in% c("continuous", "count"))
  dfq <- function(q, code) if (code == "UN") q + dfIW else df0     # IW: df = q + dfIW > q - 1 (proper)
  if (is.null(covTgxe)) covTgxe <- covT
  dfA <- dfq(Q, covT); dfT <- dfq(Q, covTgxe); dfE <- dfq(IE, covE); dfR <- nfr + dfIW
  if (is.null(eigG)) eigG <- eigen(G, symmetric = TRUE)
  if (is.null(DT)) DT <- as.matrix(dist(seq_len(Q)))
  if (is.null(DE)) DE <- as.matrix(dist(seq_len(I)))
  tc <- if (is.null(traceCells)) integer(0) else as.integer(traceCells) - 1L
  # proper prior of the hyper-mean of the intercepts: centred at the observed mean of each latent
  # dimension (0 for standardized continuous traits and for liabilities) with a very large variance
  if (is.null(m0)) m0 <- unlist(lapply(seq_len(T), function(t) {
    o <- is.finite(Ym[, t]); if (types[t] == "count") rep(mean(log(Ym[o, t] + 0.5) - Om[o, t]), nd[t]) else rep(0, nd[t]) }))
  if (is.null(vm0)) vm0 <- 1e4 * (s0A + s0T) / 0.5
  fit <- bmtmemix_gibbs(Ym, tcode, as.integer(Kcat), Om, eigG$vectors, pmax(eigG$values, 1e-8),
                        J, I, as.integer(gxe == "main"), as.integer(COV_CODES[covT]), as.integer(COV_CODES[covE]),
                        as.integer(nfT), as.integer(nfE), DT, DE, s0A, s0T, s0R,
                        dfA, dfT, dfE, dfR, df0, m0, vm0, gammaMax, initDisp,
                        as.integer(nIter), as.integer(burnIn), as.integer(thin), as.integer(nPredDraws),
                        as.integer(tc), verbose, pxScale,
                        if (priorMatch) dfIW + 3 else -1, if (priorMatch) dfIW + 3 else -1,
                        !is.null(fixCov),
                        if (is.null(fixCov)) matrix(0, 1, 1) else as.matrix(fixCov$SigmaA),
                        if (is.null(fixCov)) matrix(0, 1, 1) else as.matrix(fixCov$SigmaT),
                        if (is.null(fixCov)) matrix(0, 1, 1) else as.matrix(fixCov$Omega),
                        if (is.null(fixCov)) matrix(0, 1, 1) else as.matrix(fixCov$R),
                        if (is.null(drawRecords)) integer(0) else as.integer(cell[drawRecords] - 1L),
                        if (covTgxe == covT) -1L else as.integer(COV_CODES[covTgxe]))
  bt <- function(P) { for (t in seq_len(T)) if (types[t] == "continuous") P[, t] <- P[, t] * csd[t] + cmean[t]; P }
  mkp <- function(A) lapply(seq_len(T), function(t)
    if (types[t] %in% c("nominal", "ordinal")) { P <- matrix(A[, seq_len(Kcat[t]), t], ncol = Kcat[t]); colnames(P) <- levs[[t]]; P } else NULL)
  latnames <- unlist(lapply(seq_len(T), function(t) if (nd[t] == 1) names(Y)[t] else paste0(names(Y)[t], "_u", seq_len(nd[t]))))
  for (nm in c("SigmaA", "SigmaT", "SigmaG", "R")) dimnames(fit[[nm]]) <- list(latnames, latnames)
  dimnames(fit$SigmaE) <- list(envs, envs)
  st <- cumsum(c(1, nd))[seq_len(T)]
  gv <- fit$g[, st, drop = FALSE]; colnames(gv) <- names(Y)
  for (t in seq_len(T)) if (types[t] == "continuous") gv[, t] <- gv[, t] * csd[t]
  # predictive standard deviation of the phenotype (continuous and count traits; posterior predictive, conditional)
  psd <- sqrt(pmax(fit$pred2 - fit$pred^2, 0)); for (t in seq_len(T)) { if (types[t] == "continuous") psd[, t] <- psd[, t] * csd[t]; if (!types[t] %in% c("continuous", "count")) psd[, t] <- NA }
  structure(list(pred = bt(fit$pred), pred_sys = bt(fit$pred_sys), pred_sd = psd, gv = gv, probs = mkp(fit$prob), probs_sys = mkp(fit$prob_sys),
                 eta = fit$eta, g = fit$g, SigmaA = fit$SigmaA, SigmaT = fit$SigmaT, SigmaG = fit$SigmaG,
                 SigmaE = fit$SigmaE, Omega = fit$Omega, Ktot = fit$Ktot, R = fit$R, mu = fit$mu,
                 thresholds = lapply(seq_len(T), function(t) if (types[t] == "ordinal") fit$thresholds[t, seq_len(Kcat[t] + 1)][-c(1, Kcat[t] + 1)] else NULL),
                 traces = fit$traces, main = isTRUE(fit$useMain), gxe = gxe, types = types, Kcat = Kcat, levels = levs,
                 lines = lines, envs = envs, cell = cell, latnames = latnames, covT = covT, covTgxe = covTgxe, covE = covE, nd = nd,
                 nSave = fit$nSave, cmean = cmean, csd = csd, offsetMat = Om,
                 draws = if (is.null(drawRecords)) NULL else c(fit$draws, list(records = drawRecords)),
                 priors = list(dfA = dfA, dfT = dfT, dfE = dfE, dfR = dfR, df0 = df0, R2 = R2, gammaMax = gammaMax, priorMatch = priorMatch),
                 pxScale = pxScale, fixCov = !is.null(fixCov)),
            class = "BMTMEMix")
}

print.BMTMEMix <- function(x, ...) {
  cat("BMTMEMix fit:", length(x$types), "traits (", paste(x$types, collapse = ", "), "); G x E:", x$gxe, "\n")
  cat("Structures: traits =", x$covT, "; environments =", x$covE, "; saved draws =", x$nSave, "\n")
  if (isTRUE(x$main)) { cat("Main-effect genetic trait correlations (Sigma_A):\n"); print(round(cov2cor(x$SigmaA), 3)) }
  cat("Average within-environment genetic trait correlations (Sigma_G, identified):\n"); print(round(cov2cor(x$SigmaG), 3))
  invisible(x)
}

# =============================================================================
#  Responses: relative classes within each environment, thresholds from training records
# =============================================================================
#  Binary = above the within-environment median; categorical = within-environment tertiles.
#  In cross-validation the thresholds of an environment are computed from its training records
#  (TLTE, ULTE); for an environment without training records (TLUE, ULUE, CHRONO) they are computed
#  from all its records, i.e. the target is the relative position within the new environment.
make_responses <- function(D, train = NULL) {
  raw <- D$raw; s <- D$spec; env <- D$env
  tr <- if (is.null(train)) rep(TRUE, nrow(raw)) else seq_len(nrow(raw)) %in% train
  cls <- function(y, probs) {
    out <- integer(length(y))
    for (e in unique(env)) {
      i <- env == e; it <- i & tr; src <- if (sum(it) >= 10) y[it] else y[i]
      out[i] <- findInterval(y[i], quantile(src, probs, na.rm = TRUE, names = FALSE), left.open = TRUE) + 1L
    }
    out
  }
  data.frame(setNames(list(raw[[s[["cont"]]]], round(raw[[s[["count"]]]]),
                           cls(raw[[s[["cat"]]]], c(1/3, 2/3)), cls(raw[[s[["bin"]]]], 0.5) - 1L), unname(s)))
}

# =============================================================================
#  Accuracy metrics
# =============================================================================
metric_trait <- function(type, obs, pred, prob = NULL) {
  ok <- is.finite(obs); obs <- obs[ok]; pred <- pred[ok]; if (!is.null(prob)) prob <- prob[ok, , drop = FALSE]
  if (length(obs) < 3) return(c(M1 = NA, M2 = NA, M3 = NA, MAE = NA))
  if (type %in% c("continuous", "count")) {
    c(M1 = suppressWarnings(cor(obs, pred)), M2 = sqrt(mean((obs - pred)^2)) / mean(obs), M3 = NA, MAE = mean(abs(obs - pred)))
  } else if (type == "binary") {
    c(M1 = mean((pred > 0.5) == (obs == 1)), M2 = mean((pred - obs)^2), M3 = NA, MAE = NA)
  } else {
    lv <- as.numeric(colnames(prob)); K <- ncol(prob)
    cls <- lv[max.col(prob, ties.method = "first")]
    Yi <- outer(obs, lv, "==") * 1
    rps <- if (K > 1) mean(rowSums((t(apply(prob, 1, cumsum)) - t(apply(Yi, 1, cumsum)))[, -K, drop = FALSE]^2) / (K - 1)) else NA
    c(M1 = mean(cls == obs), M2 = mean(rowSums((prob - Yi)^2)), M3 = rps, MAE = NA)
  }
}

# =============================================================================
#  Cross-validation partitions
# =============================================================================
#  TLTE (CV2) random records; ULTE (CV1) line groups; TLUE (CV0) leave one environment out;
#  ULUE (CV00) every line group crossed with every environment (5 x I folds);
#  CHRONO: each year predicted only from earlier years (rice data).
make_folds <- function(scheme, line, env, nfold = 5, seed = 1, envOrder = NULL) {
  set.seed(seed)
  ul <- unique(line); ue <- if (is.null(envOrder)) sort(unique(env)) else envOrder; n <- length(line)
  lf <- sample(rep(seq_len(nfold), length.out = length(ul))); names(lf) <- ul
  folds <- list()
  if (scheme == "TLTE") {
    cf <- sample(rep(seq_len(nfold), length.out = n))
    for (f in seq_len(nfold)) folds[[f]] <- list(test = which(cf == f), drop = integer(0), id = paste0("f", f))
  } else if (scheme == "ULTE") {
    for (f in seq_len(nfold)) folds[[f]] <- list(test = which(lf[line] == f), drop = integer(0), id = paste0("f", f))
  } else if (scheme == "TLUE") {
    for (e in seq_along(ue)) folds[[e]] <- list(test = which(env == ue[e]), drop = integer(0), id = ue[e])
  } else if (scheme == "ULUE") {
    for (e in seq_along(ue)) for (f in seq_len(nfold)) {
      newl <- lf[line] == f
      folds[[length(folds) + 1]] <- list(test = which(env == ue[e] & newl), drop = which(env == ue[e] | newl), id = paste0(ue[e], "_f", f))
    }
  } else if (scheme == "CHRONO") {
    for (e in seq_along(ue)) if (e >= 3)
      folds[[length(folds) + 1]] <- list(test = which(env == ue[e]), drop = which(env %in% ue[e:length(ue)]), id = ue[e])
  }
  folds
}

#' Cross-validation driver. Each fit has its own seed: fit_seed + 1000 * fold + k
#' (k = 0 for the multi-trait model, k = t for the univariate model of trait t).
run_cv <- function(D, types, scheme, nfold = 5, seed = 2026, gxe = "main", covT = "UN", covE = "UN",
                   nIter = 2500, burnIn = 700, thin = 4, univariate = TRUE, multitrait = TRUE,
                   uni_traits = NULL, fit_seed = 2026, eigG = NULL, relabel = NULL, gaussian = FALSE) {
  folds <- make_folds(scheme, D$line, D$env, nfold, seed, D$envOrder)
  if (is.null(eigG)) eigG <- eigen(D$G, symmetric = TRUE)
  res <- list(); preds <- list()
  for (f in seq_along(folds)) {
    te <- folds[[f]]$test; mask <- union(te, folds[[f]]$drop)
    if (!length(te)) next
    train <- setdiff(seq_len(nrow(D$raw)), mask)
    Yd <- make_responses(D, train)
    if (!is.null(relabel)) Yd[[relabel$trait]] <- relabel$map[Yd[[relabel$trait]]]
    Ytr <- Yd; Ytr[mask, ] <- NA
    ftypes <- if (gaussian) rep("continuous", length(types)) else types
    models <- if (multitrait) list(BMTMEMix = seq_along(types)) else list()
    if (univariate) for (t in (if (is.null(uni_traits)) seq_along(types) else uni_traits)) models[[paste0("UNI_", names(Yd)[t])]] <- t
    for (m in names(models)) {
      tt <- models[[m]]; t0 <- Sys.time()
      set.seed(fit_seed + 1000 * f + (if (m == "BMTMEMix") 0 else tt))
      fit <- BMTMEMix(Ytr[, tt, drop = FALSE], ftypes[tt], D$line, D$env, D$G, gxe = gxe, covT = covT, covE = covE,
                      nIter = nIter, burnIn = burnIn, thin = thin, eigG = eigG, envLevels = D$envOrder)
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      for (k in seq_along(tt)) {
        t <- tt[k]; nm <- names(Yd)[t]
        pr <- fit$pred[fit$cell, k]
        pb <- if (types[t] %in% c("nominal", "ordinal") && !gaussian) fit$probs[[k]][fit$cell, , drop = FALSE] else NULL
        if (gaussian) {                      # Gaussian approximation: class probabilities from N(mean, residual variance)
          sdv <- sqrt(fit$R[k, k]) * fit$csd[k]
          if (types[t] == "binary") pr <- pmin(pmax(pr, 0), 1)
          if (types[t] == "count") pr <- pmax(pr, 0)
          if (types[t] %in% c("nominal", "ordinal")) {
            lv <- sort(unique(Yd[[t]])); cuts <- c(-Inf, head(lv, -1) + 0.5, Inf)
            pb <- sapply(seq_along(lv), function(c) pnorm((cuts[c + 1] - pr) / sdv) - pnorm((cuts[c] - pr) / sdv)); colnames(pb) <- lv
          }
        }
        model <- if (m == "BMTMEMix") (if (gaussian) "Gaussian" else "BMTMEMix") else "Univariate"
        for (e in unique(D$env[te])) {
          idx <- te[D$env[te] == e]
          mt <- metric_trait(types[t], Yd[idx, t], pr[idx], if (is.null(pb)) NULL else pb[idx, , drop = FALSE])
          res[[length(res) + 1]] <- data.frame(scheme = scheme, fold = f, fold_id = folds[[f]]$id, env = e, trait = nm, type = types[t],
                                               model = model, M1 = mt[["M1"]], M2 = mt[["M2"]], M3 = mt[["M3"]], MAE = mt[["MAE"]],
                                               n = length(idx), secs = secs)
        }
        pdf <- data.frame(scheme = scheme, fold = f, model = model, trait = nm, type = types[t], line = D$line[te], env = D$env[te],
                          obs = Yd[te, t], pred = pr[te], pred_sys = fit$pred_sys[fit$cell, k][te],
                          pred_sd = if (gaussian) NA else fit$pred_sd[fit$cell, k][te], p1 = NA, p2 = NA, p3 = NA)
        if (!is.null(pb)) for (c in seq_len(min(3, ncol(pb)))) pdf[[paste0("p", c)]] <- pb[te, c]
        preds[[length(preds) + 1]] <- pdf
      }
    }
  }
  list(metrics = do.call(rbind, res), preds = do.call(rbind, preds))
}
