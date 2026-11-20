# =============================================================================
#  v4_helpers.R — functions shared by the version-4 analyses (R/11–R/19, R/20)
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

V4DIR <- "outputs/v4"
dir.create(V4DIR, showWarnings = FALSE, recursive = TRUE)

# ---- draw-level quantities of a fit (fit$draws; records in fit$draws$records) ----------------------------
# latent internal scale -> trait scale for continuous traits
.bt_mean <- function(fit, t, m) if (fit$types[t] == "continuous") m * fit$csd[t] + fit$cmean[t] else m
.bt_sd   <- function(fit, t, s) if (fit$types[t] == "continuous") s * fit$csd[t] else s

# Monte Carlo standard error and posterior SD of the systematic linear predictor (latent scale) of every saved record
mc_eta <- function(fit, t) {
  nd <- fit$draws$nDraw; X <- fit$draws$sysMean[, t, seq_len(nd), drop = FALSE][, 1, ]
  if (is.null(dim(X))) X <- matrix(X, nrow = 1)
  ess <- apply(X, 1, function(x) tryCatch(posterior::ess_mean(x), error = function(e) NA_real_))
  data.table(sd_eta = apply(X, 1, sd), ess_eta = ess, mcse_eta = apply(X, 1, sd) / sqrt(pmax(ess, 1)))
}

# predictive interval (quantiles) of the posterior predictive mixture, continuous and count traits
#   cond = TRUE uses the conditional moments (phenotypic predictor), FALSE the systematic ones
pred_quantiles <- function(fit, t, probs = c(0.05, 0.95), cond = TRUE) {
  nd <- fit$draws$nDraw
  M <- (if (cond) fit$draws$condMean else fit$draws$sysMean)[, t, seq_len(nd), drop = FALSE][, 1, ]
  S <- (if (cond) fit$draws$condSD else fit$draws$sysSD)[, t, seq_len(nd), drop = FALSE][, 1, ]
  if (is.null(dim(M))) { M <- matrix(M, nrow = 1); S <- matrix(S, nrow = 1) }
  out <- matrix(NA_real_, nrow(M), length(probs))
  if (fit$types[t] == "continuous") {
    M <- .bt_mean(fit, t, M); S <- .bt_sd(fit, t, S)
    for (r in seq_len(nrow(M))) {
      Fm <- function(x) mean(pnorm((x - M[r, ]) / S[r, ]))
      lo <- min(M[r, ] - 6 * S[r, ]); hi <- max(M[r, ] + 6 * S[r, ])
      out[r, ] <- sapply(probs, function(p) uniroot(function(x) Fm(x) - p, c(lo, hi), tol = 1e-6)$root)
    }
  } else if (fit$types[t] == "count") {
    rec <- fit$draws$records; o <- fit$offsetMat[fit$cell[rec], t]
    for (r in seq_len(nrow(M))) {
      ymax <- max(5, ceiling(exp(max(M[r, ] + 6 * S[r, ]) + o[r])))
      yy <- 0:ymax
      Fy <- sapply(yy, function(y) mean(pnorm((log(y + 1) - o[r] - M[r, ]) / S[r, ])))
      out[r, ] <- sapply(probs, function(p) yy[which(Fy >= p)[1]])
    }
  }
  colnames(out) <- paste0("q", 100 * probs); out
}

# class probabilities of a Gaussian (continuous-declared) trait computed draw by draw:
# P(Y = c) = mean over draws of Phi((cut_c - m_d)/s_d) - Phi((cut_{c-1} - m_d)/s_d), cuts -Inf, 1.5, 2.5, Inf
gauss_probs_draws <- function(fit, t, levels, cond = TRUE) {
  nd <- fit$draws$nDraw
  M <- (if (cond) fit$draws$condMean else fit$draws$sysMean)[, t, seq_len(nd), drop = FALSE][, 1, ]
  S <- (if (cond) fit$draws$condSD else fit$draws$sysSD)[, t, seq_len(nd), drop = FALSE][, 1, ]
  if (is.null(dim(M))) { M <- matrix(M, nrow = 1); S <- matrix(S, nrow = 1) }
  M <- .bt_mean(fit, t, M); S <- .bt_sd(fit, t, S)
  cuts <- c(-Inf, head(levels, -1) + 0.5, Inf)
  P <- sapply(seq_along(levels), function(c) rowMeans(pnorm((cuts[c + 1] - M) / S) - pnorm((cuts[c] - M) / S)))
  if (is.null(dim(P))) P <- matrix(P, nrow = 1)
  colnames(P) <- levels; P
}

# ---- metrics on pooled out-of-fold predictions: within environment, then averaged over environments -----
mvals <- function(o, p, pr, role) {
  if (role %in% c("continuous", "count")) c(suppressWarnings(cor(o, p)), sqrt(mean((o - p)^2)) / mean(o), NA_real_)
  else if (role == "binary") c(mean((p > 0.5) == (o == 1)), mean((p - o)^2), NA_real_)
  else { Yi <- cbind(o == 1, o == 2, o == 3) * 1
    c(mean(max.col(pr, ties.method = "first") == o), mean(rowSums((pr - Yi)^2)),
      mean(((pr[, 1] - Yi[, 1])^2 + (pr[, 1] + pr[, 2] - Yi[, 1] - Yi[, 2])^2) / 2)) }
}
role_of <- function(type) ifelse(type %in% c("nominal", "ordinal"), "categorical", type)
envmetrics2 <- function(w, rows = NULL) {
  if (!is.null(rows)) w <- w[rows]
  rl <- w$role[1]
  v <- sapply(split(seq_len(nrow(w)), w$env), function(ix) {
    x <- w[ix]
    c(mvals(x$obs, x$A, cbind(x$A_p1, x$A_p2, x$A_p3), rl), mvals(x$obs, x$B, cbind(x$B_p1, x$B_p2, x$B_p3), rl)) })
  rowMeans(matrix(v, nrow = 6), na.rm = TRUE)
}
# paired cluster bootstrap over lines of B - A (w: one row per record with obs, role, env, line, A*, B*)
boot_pair <- function(w, B = as.integer(Sys.getenv("BMTME_BOOT", "1000")), seed = 2026) {
  w <- as.data.table(w); set.seed(seed)
  lines <- unique(w$line); idx <- split(seq_len(nrow(w)), w$line)
  e0 <- envmetrics2(w)
  bs <- replicate(B, { e <- envmetrics2(w, unlist(idx[sample(lines, replace = TRUE)], use.names = FALSE)); e[4:6] - e[1:3] })
  data.table(metric = c("M1", "M2", "M3"), A = e0[1:3], B = e0[4:6], diff = e0[4:6] - e0[1:3],
             lo = apply(bs, 1, quantile, 0.025, na.rm = TRUE), hi = apply(bs, 1, quantile, 0.975, na.rm = TRUE), nrec = nrow(w))
}
# wide table of two models (same records) from a long prediction table
widen2 <- function(dt, A, B, keys = c("fold", "line", "env")) {
  a <- dt[model == A, c(keys, "obs", "role", "pred", "p1", "p2", "p3"), with = FALSE]; setnames(a, c("pred", "p1", "p2", "p3"), c("A", "A_p1", "A_p2", "A_p3"))
  b <- dt[model == B, c(keys, "pred", "p1", "p2", "p3"), with = FALSE]; setnames(b, c("pred", "p1", "p2", "p3"), c("B", "B_p1", "B_p2", "B_p3"))
  merge(a, b, by = keys)
}
