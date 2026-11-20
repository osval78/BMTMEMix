# =============================================================================
# 02c_sbc_extra.R — additional simulation-based calibration (v4)
# Usage: Rscript R/02c_sbc_extra.R <nominal|diag> <first_rep> <last_rep> [px]
#  nominal: J = 40, I = 3, three traits — continuous, three-class NOMINAL (two utilities), count (Q = 4) —
#           unstructured trait matrices (the nominal model was not covered by the SBC of version 3);
#  diag   : the four traits of R/02b_sbc.R (continuous, binary, three-class ordinal, count) with DIAGONAL trait
#           matrices Sigma_A and Sigma_T (a non-UN structure; every variance drawn from its scaled inverse chi-square prior).
#  px: run with the scale move of Step 6 (only for unstructured trait matrices).
#  Parameters are drawn from the proper priors with explicit scales, the data are simulated, the model is fitted,
#  and the rank of each true value among 99 thinned draws is recorded.
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R"); source("R/v4_helpers.R")
args <- commandArgs(TRUE); cfgn <- args[1]; r1 <- as.integer(args[2]); r2 <- as.integer(args[3]); px <- length(args) > 3 && args[4] == "px"
d0 <- readRDS("data/prepared/EYT_4.rds")
J <- 40; I <- 3; Q <- 4; IE <- I - 1
lines <- rownames(d0$G)[1:J]; G <- d0$G[lines, lines]; eg <- eigen(G, symmetric = TRUE); LG <- t(chol(G))
types <- if (cfgn == "nominal") c("continuous", "nominal", "count") else c("continuous", "binary", "ordinal", "count")
covT <- if (cfgn == "diag") "DIAG" else "UN"
ps <- list(s0A = c(0.3, 0.3, 0.3, 0.05), s0T = c(0.2, 0.2, 0.2, 0.03), s0R = c(0.5, 1, 1, 0.1))
m0 <- c(0, 0, 0.5, log(8)); if (cfgn == "nominal") m0 <- c(0, 0.2, 0, log(8)); vm0 <- c(0.25, 0.25, 0.25, 0.05)
dfIW <- 4; df0 <- 5; gammaMax <- 4
riw <- function(df, S) solve(rWishart(1, df, solve(S))[, , 1])
rsic <- function(df, s) df * s / rchisq(1, df)
Cb <- local({ C <- matrix(0, I, IE); for (k in 1:IE) { nk <- sqrt(k * (k + 1)); C[1:k, k] <- 1 / nk; C[k + 1, k] <- -k / nk }; C })
out <- NULL
for (rep in r1:r2) {
  set.seed(7000 + rep + 1000 * (cfgn == "diag"))
  repeat {
    if (covT == "UN") { SA <- riw(Q + dfIW, (2 * Q + dfIW + 1) * diag(ps$s0A)); ST <- riw(Q + dfIW, (2 * Q + dfIW + 1) * diag(ps$s0T)) }
    else { SA <- diag(sapply(ps$s0A, function(s) rsic(df0, s))); ST <- diag(sapply(ps$s0T, function(s) rsic(df0, s))) }
    Om <- riw(IE + dfIW, (2 * IE + dfIW + 1) * diag(IE))
    nf <- 2; Psi <- riw(nf + dfIW, (2 * nf + dfIW + 1) * diag(ps$s0R[c(1, 4)]))
    Bm <- t(chol(Psi)) %*% matrix(rnorm(nf * 2), nf, 2) * sqrt(10)
    R <- diag(Q); fi <- c(1, 4); ci <- c(2, 3)
    R[fi, fi] <- Psi + Bm %*% t(Bm); R[fi, ci] <- Bm; R[ci, fi] <- t(Bm)
    s0Mu <- (ps$s0A + ps$s0T) / 0.5
    m <- rnorm(Q, m0, sqrt(vm0)); tau <- df0 * s0Mu / rchisq(Q, df0)
    mu <- t(sapply(1:I, function(i) rnorm(Q, m, sqrt(tau))))
    gam2 <- runif(1, 0, gammaMax)
    A <- LG %*% matrix(rnorm(J * Q), J, Q) %*% chol(SA)
    W <- array(0, c(J, IE, Q)); LO <- t(chol(Om))
    for (s in 1:Q) W[, , s] <- LG %*% matrix(rnorm(J * IE), J, IE) %*% t(LO)
    W <- array(matrix(W, J * IE, Q) %*% chol(ST), c(J, IE, Q))
    g <- array(0, c(J, I, Q)); for (s in 1:Q) g[, , s] <- matrix(A[, s], J, I) + W[, , s] %*% t(Cb)
    L <- array(0, c(J, I, Q)); for (s in 1:Q) L[, , s] <- matrix(mu[, s], J, I, byrow = TRUE) + g[, , s]
    L <- matrix(L, J * I, Q) + matrix(rnorm(J * I * Q), J * I, Q) %*% chol(R)
    if (cfgn == "nominal") {
      cl <- ifelse(pmax(L[, 2], L[, 3]) < 0, 3L, ifelse(L[, 2] > L[, 3], 1L, 2L))
      Yd <- data.frame(y1 = L[, 1], y2 = cl, y4 = floor(exp(L[, 4])))
      ok <- all(table(factor(cl, 1:3)) >= 3) && max(Yd$y4) < 1e6
    } else {
      Yd <- data.frame(y1 = L[, 1], y2 = as.integer(L[, 2] > 0), y3 = findInterval(L[, 3], c(0, gam2)) + 1L, y4 = floor(exp(L[, 4])))
      ok <- length(unique(Yd$y3)) == 3 && length(unique(Yd$y2)) == 2 && max(Yd$y4) < 1e6
    }
    if (ok) break
  }
  line <- rep(lines, I); env <- rep(paste0("E", 1:I), each = J)
  sd1 <- sd(Yd$y1); mn1 <- mean(Yd$y1)
  pr <- ps; pr$s0A[1] <- ps$s0A[1] / sd1^2; pr$s0T[1] <- ps$s0T[1] / sd1^2; pr$s0R[1] <- ps$s0R[1] / sd1^2
  m0s <- m0; m0s[1] <- (m0[1] - mn1) / sd1; vm0s <- vm0; vm0s[1] <- vm0[1] / sd1^2
  set.seed(rep)
  fit <- BMTMEMix(Yd, types, line, env, G, covT = covT, nIter = 3000, burnIn = 1020, thin = 20, eigG = eg, priorScales = pr,
                  m0 = m0s, vm0 = vm0s, dfIW = dfIW, df0 = df0, gammaMax = gammaMax, traceCells = 1, pxScale = px)
  tr <- fit$traces; ix <- function(r, c) (c - 1) * Q + r
  cE <- sum(diag(Om)) / IE
  if (cfgn == "nominal") {
    truth <- c(SigmaA11 = SA[1, 1] / sd1^2, SigmaA22_nom = SA[2, 2], corA23_nom = cov2cor(SA)[2, 3], corA12 = cov2cor(SA)[1, 2],
               SigmaTn22_nom = ST[2, 2] * cE, R11 = R[1, 1] / sd1^2, R44 = R[4, 4], mu_nom1 = mu[1, 2], eta_cell1_nom1 = mu[1, 2] + g[1, 1, 2])
    draws <- cbind(SigmaA11 = tr$SigmaA[, ix(1, 1)], SigmaA22_nom = tr$SigmaA[, ix(2, 2)],
                   corA23_nom = tr$SigmaA[, ix(2, 3)] / sqrt(tr$SigmaA[, ix(2, 2)] * tr$SigmaA[, ix(3, 3)]),
                   corA12 = tr$SigmaA[, ix(1, 2)] / sqrt(tr$SigmaA[, ix(1, 1)] * tr$SigmaA[, ix(2, 2)]),
                   SigmaTn22_nom = tr$SigmaT[, ix(2, 2)], R11 = tr$R[, ix(1, 1)], R44 = tr$R[, ix(4, 4)],
                   mu_nom1 = tr$mu[, 1 * I + 1], eta_cell1_nom1 = tr$eta[, 2])
  } else {
    k <- which(colSums(abs(tr$thresholds)) > 0)
    truth <- c(SigmaA11 = SA[1, 1] / sd1^2, SigmaA22_bin = SA[2, 2], SigmaA44 = SA[4, 4], SigmaTn33 = ST[3, 3] * cE,
               R11 = R[1, 1] / sd1^2, R44 = R[4, 4], gamma2 = gam2, mu_count = mu[1, 4], eta_cell1_bin = mu[1, 2] + g[1, 1, 2])
    draws <- cbind(SigmaA11 = tr$SigmaA[, ix(1, 1)], SigmaA22_bin = tr$SigmaA[, ix(2, 2)], SigmaA44 = tr$SigmaA[, ix(4, 4)],
                   SigmaTn33 = tr$SigmaT[, ix(3, 3)], R11 = tr$R[, ix(1, 1)], R44 = tr$R[, ix(4, 4)], gamma2 = tr$thresholds[, k[1]],
                   mu_count = tr$mu[, 3 * I + 1], eta_cell1_bin = tr$eta[, 2])
  }
  draws <- draws[seq_len(min(99, nrow(draws))), , drop = FALSE]
  out <- rbind(out, data.frame(config = cfgn, px = px, rep = rep, parameter = names(truth), rank = colSums(sweep(draws, 2, truth, "<")), L = nrow(draws)))
  if (rep %% 20 == 0) cat("sbc", rep, "\n")
}
dir.create(file.path(V4DIR, "sbc"), showWarnings = FALSE, recursive = TRUE)
fwrite(out, file.path(V4DIR, "sbc", sprintf("sbc_%s%s_%03d_%03d.csv", cfgn, if (px) "px" else "", r1, r2)))
cat("done\n")
