# =============================================================================
# 02b_sbc.R — simulation-based calibration (Talts et al., 2018) of the sampler on a small problem
# Usage: Rscript R/02b_sbc.R <first_rep> <last_rep> [px]
#   px: run the sampler with the scale move of Step 6 (pxScale = TRUE); output file sbcpx_*.csv
# Parameters drawn from the (proper) priors of BMTMEMix with explicit prior scales; data with J = 40 lines,
# I = 3 environments and four traits (continuous, binary, ordinal with K = 3, count); main + contrast G x E.
# For each replicate the rank of the true value among L = 99 thinned posterior draws is recorded; under a
# correct sampler the ranks are uniform on {0, ..., 99}.
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R")
args <- commandArgs(TRUE); r1 <- as.integer(args[1]); r2 <- as.integer(args[2]); px <- length(args) > 2 && args[3] == "px"
d0 <- readRDS("data/prepared/EYT_4.rds")
J <- 40; I <- 3; Q <- 4; IE <- I - 1
lines <- rownames(d0$G)[1:J]; G <- d0$G[lines, lines]; eg <- eigen(G, symmetric = TRUE)
types <- c("continuous", "binary", "ordinal", "count")
ps <- list(s0A = c(0.3, 0.3, 0.3, 0.05), s0T = c(0.2, 0.2, 0.2, 0.03), s0R = c(0.5, 1, 1, 0.1))
m0 <- c(0, 0, 0.5, log(8)); vm0 <- c(0.25, 0.25, 0.25, 0.05)
dfIW <- 4; df0 <- 5; gammaMax <- 4
riw <- function(df, S) solve(rWishart(1, df, solve(S))[, , 1])
Cb <- local({ C <- matrix(0, I, IE); for (k in 1:IE) { nk <- sqrt(k * (k + 1)); C[1:k, k] <- 1 / nk; C[k + 1, k] <- -k / nk }; C })
out <- NULL
for (rep in r1:r2) {
  set.seed(5000 + rep)
  repeat {
    SA <- riw(Q + dfIW, (2 * Q + dfIW + 1) * diag(ps$s0A)); ST <- riw(Q + dfIW, (2 * Q + dfIW + 1) * diag(ps$s0T))
    Om <- riw(IE + dfIW, (2 * IE + dfIW + 1) * diag(IE))
    nf <- 2; Psi <- riw(nf + dfIW, (2 * nf + dfIW + 1) * diag(ps$s0R[c(1, 4)]))
    Bm <- t(chol(Psi)) %*% matrix(rnorm(nf * 2), nf, 2) * sqrt(10)          # B | Psi ~ MN(0, Psi, 10 I)
    R <- diag(Q); fi <- c(1, 4); ci <- c(2, 3)
    R[fi, fi] <- Psi + Bm %*% t(Bm); R[fi, ci] <- Bm; R[ci, fi] <- t(Bm)
    s0Mu <- (ps$s0A + ps$s0T) / 0.5
    m <- rnorm(Q, m0, sqrt(vm0)); tau <- df0 * s0Mu / rchisq(Q, df0)
    mu <- t(sapply(1:I, function(i) rnorm(Q, m, sqrt(tau))))
    gam2 <- runif(1, 0, gammaMax)
    LG <- t(chol(G))
    A <- LG %*% matrix(rnorm(J * Q), J, Q) %*% chol(SA)
    W <- array(0, c(J, IE, Q)); LO <- t(chol(Om))
    for (s in 1:Q) W[, , s] <- LG %*% matrix(rnorm(J * IE), J, IE) %*% t(LO)
    W <- array(matrix(W, J * IE, Q) %*% chol(ST), c(J, IE, Q))
    g <- array(0, c(J, I, Q)); for (s in 1:Q) g[, , s] <- matrix(A[, s], J, I) + W[, , s] %*% t(Cb)
    L <- array(0, c(J, I, Q)); for (s in 1:Q) L[, , s] <- matrix(mu[, s], J, I, byrow = TRUE) + g[, , s]
    L <- matrix(L, J * I, Q) + matrix(rnorm(J * I * Q), J * I, Q) %*% chol(R)
    Yd <- data.frame(y1 = L[, 1], y2 = as.integer(L[, 2] > 0), y3 = findInterval(L[, 3], c(0, gam2)) + 1L, y4 = floor(exp(L[, 4])))
    if (length(unique(Yd$y3)) == 3 && length(unique(Yd$y2)) == 2 && max(Yd$y4) < 1e6) break
  }
  line <- rep(lines, I); env <- rep(paste0("E", 1:I), each = J)
  # the continuous trait is standardized internally: express its truths on that scale
  sd1 <- sd(Yd$y1); mn1 <- mean(Yd$y1)
  pr <- ps; pr$s0A[1] <- ps$s0A[1] / sd1^2; pr$s0T[1] <- ps$s0T[1] / sd1^2; pr$s0R[1] <- ps$s0R[1] / sd1^2
  m0s <- m0; m0s[1] <- (m0[1] - mn1) / sd1; vm0s <- vm0; vm0s[1] <- vm0[1] / sd1^2
  set.seed(rep)
  fit <- BMTMEMix(Yd, types, line, env, G, nIter = 3000, burnIn = 1020, thin = 20, eigG = eg, priorScales = pr,
                  m0 = m0s, vm0 = vm0s, dfIW = dfIW, df0 = df0, gammaMax = gammaMax, traceCells = 1, pxScale = px)
  tr <- fit$traces; k <- which(colSums(abs(tr$thresholds)) > 0)
  cE <- sum(diag(Om)) / IE
  truth <- c(SigmaA11 = SA[1, 1] / sd1^2, SigmaA44 = SA[4, 4], corA14 = cov2cor(SA)[1, 4], SigmaTn33 = ST[3, 3] * cE,
             R11 = R[1, 1] / sd1^2, R44 = R[4, 4], gamma2 = gam2, mu_count = mu[1, 4],
             eta_cell1_bin = mu[1, 2] + g[1, 1, 2])
  draws <- cbind(SigmaA11 = tr$SigmaA[, 1], SigmaA44 = tr$SigmaA[, 16],
                 corA14 = tr$SigmaA[, 13] / sqrt(tr$SigmaA[, 1] * tr$SigmaA[, 16]), SigmaTn33 = tr$SigmaT[, 11],
                 R11 = tr$R[, 1], R44 = tr$R[, 16], gamma2 = tr$thresholds[, k[1]], mu_count = tr$mu[, 3 * I + 1],
                 eta_cell1_bin = tr$eta[, 2])
  draws <- draws[seq_len(min(99, nrow(draws))), , drop = FALSE]
  out <- rbind(out, data.frame(rep = rep, parameter = names(truth), rank = colSums(sweep(draws, 2, truth, "<")), L = nrow(draws)))
  if (rep %% 20 == 0) cat("sbc", rep, "\n")
}
dir.create("outputs/simulation", showWarnings = FALSE, recursive = TRUE)
write.csv(out, sprintf("outputs/simulation/%s_%03d_%03d.csv", if (px) "sbcpx" else "sbc", r1, r2), row.names = FALSE)
