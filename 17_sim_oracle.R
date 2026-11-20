# =============================================================================
# 17_sim_oracle.R — ideal (known covariances) versus estimated multi-trait and univariate models (v4)
# Usage: Rscript R/17_sim_oracle.R <config 1-12> <first_rep> <last_rep>
#  J = 300 lines (genomic relationships of the first 300 lines of EYT_4), I = 3 environments, two traits:
#  a target (continuous or three-class ordinal) and a continuous secondary trait of high signal (h2 = 0.7).
#  Configurations: target type x target signal (h2 = 0.2 or 0.6, liability scale for the ordinal target)
#  x genetic correlation rho (0, 0.5, 0.8; the same for the main and the G x E effects); residual correlation 0.3.
#  20% of the lines are untested; two masks: coincident (all traits of the test lines masked) and target-specific
#  (only the target masked; the secondary trait of the test lines observed in every environment).
#  Four fits per mask: MT-true (Sigma_A, Sigma_T, Omega and R fixed at their true values), UNI-true (univariate model
#  with the true marginal variances of the target), MT-est and UNI-est (all covariances estimated; default priors).
#  Accuracy = correlation between predicted and true genetic values (g = a + u) of the test lines within each
#  environment, averaged over environments.
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R"); source("R/v4_helpers.R")
args <- commandArgs(TRUE); cfg <- as.integer(args[1]); r1 <- as.integer(args[2]); r2 <- as.integer(args[3])
grid <- expand.grid(rho = c(0, 0.5, 0.8), h2 = c(0.2, 0.6), ttype = c("continuous", "ordinal"), stringsAsFactors = FALSE)
cf <- grid[cfg, ]
d0 <- readRDS("data/prepared/EYT_4.rds"); J <- 300; I <- 3
lines <- rownames(d0$G)[1:J]; G <- d0$G[lines, lines]; eg <- eigen(G, symmetric = TRUE); LG <- t(chol(G))
Cb <- local({ C <- matrix(0, I, I - 1); for (k in 1:(I - 1)) { nk <- sqrt(k * (k + 1)); C[1:k, k] <- 1 / nk; C[k + 1, k] <- -k / nk }; C })
Om <- matrix(c(1.1, 0.25, 0.25, 0.9), 2, 2)                 # trace / (I - 1) = 1
Cr <- matrix(c(1, cf$rho, cf$rho, 1), 2, 2)
Vg <- c(if (cf$ttype == "ordinal") cf$h2 / (1 - cf$h2) else cf$h2, 0.7)
Rv <- c(if (cf$ttype == "ordinal") 1 else 1 - cf$h2, 0.3)
SA <- 0.6 * diag(sqrt(Vg)) %*% Cr %*% diag(sqrt(Vg)); ST <- 0.6 * diag(sqrt(Vg)) %*% Cr %*% diag(sqrt(Vg))   # Sigma_G = SA + (2/3) ST
Rt <- diag(Rv); Rt[1, 2] <- Rt[2, 1] <- 0.3 * sqrt(Rv[1] * Rv[2])
types <- c(cf$ttype, "continuous")
line <- rep(lines, I); env <- rep(paste0("E", 1:I), each = J); cell <- seq_len(J * I)
RES <- NULL
for (rep in r1:r2) {
  set.seed(30000 + 100 * cfg + rep)
  A <- LG %*% matrix(rnorm(J * 2), J, 2) %*% chol(SA)
  U <- array(0, c(J, I, 2)); Wl <- LG %*% matrix(rnorm(J * 2 * (I - 1)), J, 2 * (I - 1))   # columns: (trait, contrast)
  Kc <- kronecker(ST, Om); Wm <- Wl %*% chol(Kc)                  # J x (2 * (I-1)), Sigma_T (x) Omega across columns
  for (s in 1:2) U[, , s] <- Wm[, (s - 1) * (I - 1) + 1:(I - 1)] %*% t(Cb)
  g <- sapply(1:2, function(s) as.vector(A[, s] + U[, , s]))      # (J*I) x 2, lines fastest
  mu <- cbind(rnorm(I, 0, 0.5), rnorm(I, 0, 0.5)); if (cf$ttype == "ordinal") mu[, 1] <- 0.45 + rnorm(I, 0, 0.2)
  E <- matrix(rnorm(J * I * 2), J * I, 2) %*% chol(Rt)
  L <- g + mu[rep(1:I, each = J), ] + E
  y1 <- if (cf$ttype == "ordinal") findInterval(L[, 1], c(0, 0.9)) + 1 else 10 + 2 * L[, 1]
  Y <- data.frame(tgt = y1, sec = 50 + 5 * L[, 2])
  test_lines <- sample(lines, J / 5); te <- which(line %in% test_lines)
  Yt <- Y; Yt$tgt[te] <- NA
  acc <- function(fit, k = 1) mean(sapply(1:I, function(i) { ix <- te[env[te] == paste0("E", i)]; cor(fit$g[fit$cell[ix], k], g[ix, 1]) }))
  # true covariances on the internal scale of each fit (continuous traits standardized by the SD of the observed data)
  fix <- function(csd) { Dm <- diag(1 / csd, length(csd)); list(SigmaA = Dm %*% SA[seq_along(csd), seq_along(csd), drop = FALSE] %*% Dm,
                                                              SigmaT = Dm %*% ST[seq_along(csd), seq_along(csd), drop = FALSE] %*% Dm, Omega = Om,
                                                              R = Dm %*% Rt[seq_along(csd), seq_along(csd), drop = FALSE] %*% Dm) }
  scale_y <- c(if (cf$ttype == "ordinal") 1 else 2, 5)              # y = a + b L
  mc <- list(nIter = 2000, burnIn = 500, thin = 3)
  for (mk in c("coincident", "target-specific")) {
    Ym <- Yt; if (mk == "coincident") Ym$sec[te] <- NA
    csd <- c(if (cf$ttype == "ordinal") 1 else sd(Ym$tgt, na.rm = TRUE) / scale_y[1], sd(Ym$sec, na.rm = TRUE) / scale_y[2])
    for (oracle in c(FALSE, TRUE)) {
      set.seed(rep + 77 * oracle); t0 <- Sys.time()
      fit <- BMTMEMix(Ym, types, line, env, G, nIter = mc$nIter, burnIn = mc$burnIn, thin = mc$thin, eigG = eg,
                      fixCov = if (oracle) fix(csd) else NULL)
      RES <- rbind(RES, data.frame(cfg = cfg, rho = cf$rho, h2 = cf$h2, ttype = cf$ttype, rep = rep, mask = mk,
                                   model = if (oracle) "MT-true" else "MT-est", acc = acc(fit),
                                   secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }
  csd <- if (cf$ttype == "ordinal") 1 else sd(Yt$tgt, na.rm = TRUE) / scale_y[1]
  for (oracle in c(FALSE, TRUE)) {
    set.seed(rep + 99 + 77 * oracle)
    fu <- BMTMEMix(Yt[, "tgt", drop = FALSE], types[1], line, env, G, nIter = mc$nIter, burnIn = mc$burnIn, thin = mc$thin, eigG = eg,
                   fixCov = if (oracle) fix(csd) else NULL)
    a <- acc(fu)
    RES <- rbind(RES, data.frame(cfg = cfg, rho = cf$rho, h2 = cf$h2, ttype = cf$ttype, rep = rep, mask = c("coincident", "target-specific"),
                                 model = if (oracle) "UNI-true" else "UNI-est", acc = a, secs = NA))
  }
  cat("rep", rep, "done", format(Sys.time()), "\n")
}
dir.create(file.path(V4DIR, "sim_oracle"), showWarnings = FALSE, recursive = TRUE)
fwrite(RES, file.path(V4DIR, "sim_oracle", sprintf("oracle_cfg%02d_%03d_%03d.csv", cfg, r1, r2)))
cat("done cfg", cfg, "\n")
