# =============================================================================
# 02_simulation_check.R — replicated simulation study (version 3)
# Usage: Rscript R/02_simulation_check.R <setting: base|hard> <first_rep> <last_rep>
#  base: J = 400, I = 3, five response types (Q = 6), balanced, moderate class frequencies
#  hard: rare classes, many zero counts, 40% of the cells missing at random, stronger residual correlations
# For every replicate: (1) full-data fit -> bias, RMSE and 90% interval coverage of identified parameters;
# (2) genetic accuracy (correlation with the true genetic values) for 20% untested lines whose target trait is
#     masked while their other traits are observed (trait-assisted, MT-TA), with all traits masked (MT-all),
#     and with the univariate model (UNI), for a continuous and an ordinal target.
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R")
args <- commandArgs(TRUE); setting <- args[1]; r1 <- as.integer(args[2]); r2 <- as.integer(args[3])
dir.create("outputs/simulation", showWarnings = FALSE, recursive = TRUE)
d0 <- readRDS("data/prepared/EYT_4.rds")
J <- 400; I <- 3; Q <- 6
lines <- rownames(d0$G)[1:J]; G <- d0$G[lines, lines]; eg <- eigen(G, symmetric = TRUE)
pn <- c("cont", "cnt", "nom1", "nom2", "bin", "ord")
types <- c("continuous", "count", "nominal", "binary", "ordinal")
cor2cov <- function(C, v) diag(sqrt(v)) %*% C %*% diag(sqrt(v))
CA <- matrix(c(1, .5, .3, .2, .3, .4,  .5, 1, .2, .1, .2, .3,  .3, .2, 1, .4, .2, .2,
               .2, .1, .4, 1, .1, .1,  .3, .2, .2, .1, 1, .3,  .4, .3, .2, .1, .3, 1), Q, Q)
SA <- cor2cov(CA, c(0.8, 0.6, 0.8, 0.8, 0.8, 0.8))
CT <- CA; CT[CT > 0 & CT < 1] <- CT[CT > 0 & CT < 1] * 0.8
ST <- cor2cov(CT, c(0.4, 0.3, 0.4, 0.4, 0.4, 0.4))
Om <- matrix(c(1.2, 0.3, 0.3, 0.8), 2, 2)                    # G x E covariance on the 2 environmental contrasts
Cb <- local({ C <- matrix(0, I, I - 1); for (k in 1:(I - 1)) { nk <- sqrt(k * (k + 1)); C[1:k, k] <- 1 / nk; C[k + 1, k] <- -k / nk }; C })
mk_setting <- function(setting) {
  Rt <- diag(Q); Rt[1, 1] <- 0.6; Rt[2, 2] <- 0.5; Rt[1, 2] <- Rt[2, 1] <- 0.2
  if (setting == "base") {
    Rt[1, 5] <- Rt[5, 1] <- 0.3; Rt[2, 6] <- Rt[6, 2] <- 0.2
    list(Rt = Rt, gam = c(0, 0.8, 1.6), mu0 = c(0, 1.5, 0.2, -0.2, 0, 0.6), miss = 0)
  } else {
    Rt[1, 5] <- Rt[5, 1] <- 0.45; Rt[2, 6] <- Rt[6, 2] <- 0.35
    list(Rt = Rt, gam = c(0, 1.5, 2.5), mu0 = c(0, -0.6, 0.6, 0.3, 1.0, -0.9), miss = 0.4)
  }
}
S <- mk_setting(setting)
# identified true values: Sigma_A, normalized Sigma_T x Omega, Sigma_G, residual covariance, thresholds
cE <- sum(diag(Om)) / (I - 1)
STn <- ST * cE; SEn <- Cb %*% (Om / cE) %*% t(Cb); SG <- SA + STn * sum(diag(SEn)) / I
tracepars <- function(fit, Dsc) {
  tr <- fit$traces; n <- nrow(tr$SigmaA); out <- list(); truth <- list(); fam <- list()
  scl <- diag(Dsc)
  addp <- function(name, draws, tv, family) { out[[name]] <<- draws; truth[[name]] <<- tv; fam[[name]] <<- family }
  for (s in 1:Q) {
    addp(sprintf("SigmaA[%s]", pn[s]), tr$SigmaA[, (s - 1) * Q + s], SA[s, s] * scl[s]^2, "Sigma_A variance")
    addp(sprintf("SigmaTn[%s]", pn[s]), tr$SigmaT[, (s - 1) * Q + s], STn[s, s] * scl[s]^2, "Sigma_T x Omega variance")
  }
  for (s in 1:(Q - 1)) for (l in (s + 1):Q) {
    addp(sprintf("corA[%s,%s]", pn[s], pn[l]), sapply(1:n, function(r) cov2cor(matrix(tr$SigmaA[r, ], Q, Q))[s, l]), cov2cor(SA)[s, l], "Sigma_A correlation")
    addp(sprintf("corG[%s,%s]", pn[s], pn[l]), sapply(1:n, function(r) cov2cor(matrix(tr$SigmaG[r, ], Q, Q))[s, l]), cov2cor(SG)[s, l], "Sigma_G correlation")
  }
  for (i in 1:I) for (k in i:I) addp(sprintf("SigmaE_GxE[%d,%d]", i, k), tr$SigmaE[, (k - 1) * I + i], SEn[i, k], "G x E environmental covariance")
  for (s in c(1, 2)) addp(sprintf("R[%s]", pn[s]), tr$R[, (s - 1) * Q + s], S$Rt[s, s] * scl[s]^2, "Residual variance")
  for (pp in list(c(1, 2), c(1, 5), c(2, 6))) addp(sprintf("corR[%s,%s]", pn[pp[1]], pn[pp[2]]),
      tr$R[, (pp[2] - 1) * Q + pp[1]] / sqrt(tr$R[, (pp[1] - 1) * Q + pp[1]] * tr$R[, (pp[2] - 1) * Q + pp[2]]),
      cov2cor(S$Rt)[pp[1], pp[2]], "Residual correlation")
  k <- which(colSums(abs(tr$thresholds)) > 0)
  addp("gamma_2", tr$thresholds[, k[1]], S$gam[2], "Ordinal threshold"); addp("gamma_3", tr$thresholds[, k[2]], S$gam[3], "Ordinal threshold")
  data.frame(parameter = names(out), family = unlist(fam), true = unlist(truth),
             mean = sapply(out, mean), q05 = sapply(out, quantile, 0.05), q95 = sapply(out, quantile, 0.95), row.names = NULL)
}
RES <- NULL; ACC <- NULL
for (rep in r1:r2) {
  set.seed(10000 * (setting == "hard") + rep)
  mu <- matrix(S$mu0, I, Q, byrow = TRUE) + matrix(rnorm(I * Q, 0, .3), I, Q)
  LG <- t(chol(G))
  A <- LG %*% matrix(rnorm(J * Q), J, Q) %*% chol(SA)
  W <- array(0, c(J, I - 1, Q)); LO <- t(chol(Om))
  for (s in 1:Q) W[, , s] <- LG %*% matrix(rnorm(J * (I - 1)), J, I - 1) %*% t(LO)
  W <- array(matrix(W, J * (I - 1), Q) %*% chol(ST), c(J, I - 1, Q))
  g <- array(0, c(J, I, Q)); for (s in 1:Q) g[, , s] <- matrix(A[, s], J, I) + W[, , s] %*% t(Cb)
  E <- matrix(rnorm(J * I * Q), J * I, Q) %*% chol(S$Rt)
  L <- array(0, c(J, I, Q)); for (s in 1:Q) L[, , s] <- matrix(mu[, s], J, I, byrow = TRUE) + g[, , s]
  L <- matrix(L, J * I, Q) + E
  gam <- c(-Inf, S$gam, Inf)
  Yd <- data.frame(cont = L[, 1], cnt = floor(exp(L[, 2])), nom = apply(cbind(L[, 3:4], 0), 1, which.max),
                   bin = as.integer(L[, 5] > 0), ord = findInterval(L[, 6], gam[2:4]) + 1L)
  line <- rep(lines, I); env <- rep(paste0("E", 1:I), each = J)
  if (S$miss > 0) { keep <- sort(sample(J * I, round((1 - S$miss) * J * I))); Yd <- Yd[keep, ]; line <- line[keep]; env <- env[keep]; Lk <- keep } else Lk <- seq_len(J * I)
  gt <- matrix(g, J * I, Q)[Lk, ]
  set.seed(rep); t0 <- Sys.time()
  fit <- BMTMEMix(Yd, types, line, env, G, nIter = 3000, burnIn = 1000, thin = 4, eigG = eg)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  Dsc <- diag(c(1 / fit$csd[1], rep(1, Q - 1)))
  tp <- tracepars(fit, Dsc); tp$setting <- setting; tp$rep <- rep; tp$secs <- secs
  tp$freq <- paste(sprintf("nom %s; ord %s; zeros %.2f", paste(table(Yd$nom), collapse = "/"), paste(table(Yd$ord), collapse = "/"), mean(Yd$cnt == 0)))
  RES <- rbind(RES, tp)
  gacc <- data.frame(setting = setting, rep = rep, dim = pn, cor_g_full = diag(cor(gt, fit$g[fit$cell, ])))
  # ---- genetic accuracy for untested lines ----
  testl <- sample(lines, round(0.2 * J)); te <- which(line %in% testl)
  for (tg in c(1, 5)) {                       # trait index in Yd: 1 = continuous, 5 = ordinal
    dim <- if (tg == 1) 1 else 6
    Yta <- Yd; Yta[te, tg] <- NA
    set.seed(rep + 1); f_ta <- BMTMEMix(Yta, types, line, env, G, nIter = 3000, burnIn = 1000, thin = 4, eigG = eg)
    Yu <- Yd[, tg, drop = FALSE]; Yu[te, 1] <- NA
    set.seed(rep + 2); f_u <- BMTMEMix(Yu, types[tg], line, env, G, nIter = 3000, burnIn = 1000, thin = 4, eigG = eg)
    ACC <- rbind(ACC, data.frame(setting = setting, rep = rep, target = pn[dim], model = c("MT-TA", "UNI"),
                                 cor_g = c(cor(gt[te, dim], f_ta$g[f_ta$cell[te], dim]), cor(gt[te, dim], f_u$g[f_u$cell[te], 1]))))
  }
  Yall <- Yd; Yall[te, ] <- NA
  set.seed(rep + 3); f_all <- BMTMEMix(Yall, types, line, env, G, nIter = 3000, burnIn = 1000, thin = 4, eigG = eg)
  ACC <- rbind(ACC, data.frame(setting = setting, rep = rep, target = c("cont", "ord"), model = "MT-all",
                               cor_g = c(cor(gt[te, 1], f_all$g[f_all$cell[te], 1]), cor(gt[te, 6], f_all$g[f_all$cell[te], 6]))))
  write.csv(gacc, sprintf("outputs/simulation/gv_%s_%03d.csv", setting, rep), row.names = FALSE)
  cat(setting, "rep", rep, "done", round(secs), "s\n")
}
write.csv(RES, sprintf("outputs/simulation/params_%s_%03d_%03d.csv", setting, r1, r2), row.names = FALSE)
write.csv(ACC, sprintf("outputs/simulation/accuracy_%s_%03d_%03d.csv", setting, r1, r2), row.names = FALSE)
