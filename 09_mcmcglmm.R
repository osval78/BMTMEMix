# Usage: Rscript R/09_mcmcglmm.R <dataset> [nfolds]
# Competitor: multi-response generalized linear mixed model fitted with MCMCglmm (Hadfield, 2010),
# scenario 2 (continuous Gaussian, Poisson count with log-normal overdispersion, ordinal and binary
# probit "threshold" traits), TLTE folds identical to BMTMEMix. Model: trait-specific environment means,
# genomic line effects with an unstructured trait covariance us(trait):Line (ginverse = G^-1), independent
# trait-specific line x environment effects, and an unstructured residual covariance whose block for the two
# threshold traits is fixed (fix = 3).
source("R/00_settings.R"); source("R/BMTMEMix.R")
suppressPackageStartupMessages({ library(MCMCglmm); library(Matrix) })
args <- commandArgs(TRUE); d <- args[1]; nf <- if (length(args) > 1) as.integer(args[2]) else NFOLD
D <- readRDS(sprintf("data/prepared/%s.rds", d)); types <- D$types_S2
folds <- make_folds("TLTE", D$line, D$env, NFOLD, FOLD_SEED)[seq_len(nf)]
Ginv <- as(solve(D$G), "dgCMatrix"); dimnames(Ginv) <- dimnames(D$G)
res <- NULL
for (f in seq_along(folds)) {
  te <- folds[[f]]$test; train <- setdiff(seq_len(nrow(D$raw)), te)
  Yd <- make_responses(D, train); tn <- names(Yd)
  dat <- data.frame(Line = factor(D$line, levels = rownames(D$G)), env = factor(D$env), LineEnv = factor(paste(D$line, D$env)))
  m1 <- mean(Yd[train, 1]); s1 <- sd(Yd[train, 1])
  dat$y1 <- (Yd[[1]] - m1) / s1; dat$y2 <- Yd[[2]]; dat$y3 <- Yd[[3]]; dat$y4 <- Yd[[4]] + 1L   # threshold traits coded 1..K
  dat[te, c("y1", "y2", "y3", "y4")] <- NA
  prior <- list(R = list(V = diag(4), nu = 5, fix = 3),
                G = list(G1 = list(V = diag(4) * 0.5, nu = 5), G2 = list(V = diag(4) * 0.25, nu = 5)))
  set.seed(FIT_SEED + 1000 * f); t0 <- Sys.time()
  m <- MCMCglmm(cbind(y1, y2, y3, y4) ~ trait - 1 + trait:env, random = ~ us(trait):Line + idh(trait):LineEnv,
                rcov = ~ us(trait):units, family = c("gaussian", "poisson", "threshold", "threshold"),
                ginverse = list(Line = Ginv), data = dat, prior = prior, pl = TRUE, verbose = FALSE,
                nitt = 2500, burnin = 500, thin = 4)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  n <- nrow(dat); Liab <- m$Liab          # draws x (n * 4), trait-major
  CP <- m$CP                              # cutpoints (draws x ...), ordered by threshold trait
  ncp3 <- length(unique(Yd[[3]])) - 2     # free cut-points of the ordinal trait (first fixed at 0)
  for (t in 1:4) {
    L <- Liab[, (t - 1) * n + te, drop = FALSE]
    if (t == 1) { pr <- colMeans(L) * s1 + m1; pb <- NULL }
    if (t == 2) { pr <- colMeans(exp(L)); pb <- NULL }
    if (t == 3) {
      cps <- if (ncp3 > 0) as.matrix(CP)[, seq_len(ncp3), drop = FALSE] else matrix(numeric(0), nrow(L), 0)
      K <- ncp3 + 2
      pb <- sapply(seq_len(K), function(c) {
        lo <- if (c == 1) -Inf else if (c == 2) 0 else cps[, c - 2]
        hi <- if (c == 1) 0 else if (c == K) Inf else cps[, c - 1]
        colMeans(L > lo & L <= hi) }); colnames(pb) <- seq_len(K); pr <- NA
    }
    if (t == 4) { pr <- colMeans(L > 0); pb <- NULL }
    for (e in unique(D$env[te])) {
      w <- which(D$env[te] == e); idx <- te[w]
      mt <- metric_trait(types[t], Yd[idx, t], if (length(pr) > 1) pr[w] else rep(NA, length(w)), if (is.null(pb)) NULL else pb[w, , drop = FALSE])
      res <- rbind(res, data.frame(dataset = d, scenario = "S2", scheme = "TLTE", fold = f, env = e, trait = tn[t], type = types[t],
                                   model = "MCMCglmm", M1 = mt[["M1"]], M2 = mt[["M2"]], M3 = mt[["M3"]], MAE = mt[["MAE"]], n = length(idx), secs = secs))
    }
  }
  cat("fold", f, round(secs), "s\n")
}
dir.create("outputs/cv/S2_mcmcglmm", showWarnings = FALSE, recursive = TRUE)
write.csv(res, sprintf("outputs/cv/S2_mcmcglmm/%s_TLTE_metrics.csv", d), row.names = FALSE)
