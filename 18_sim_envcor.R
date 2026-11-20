# =============================================================================
# 18_sim_envcor.R — untested environment whose target–auxiliary genetic correlation differs (v4)
# Usage: Rscript R/18_sim_envcor.R <first_rep> <last_rep>
#  J = 300 lines, I = 4 environments, two continuous traits (target h2 = 0.4, auxiliary h2 = 0.7, residual
#  correlation 0.3). Genetic values g_ti = sqrt(0.6) z_t + sqrt(0.4) d_ti with z ~ N(0, Cm (x) G), corr 0.3, and
#  environment-specific deviations d_.i ~ N(0, Cd_i (x) G):
#   "constant": corr(d) = 0.9 in every environment (within-environment genetic correlation 0.54 everywhere);
#   "changing": corr(d) = 0.9 in environments 1-3 and -0.9 in environment 4 (within-environment correlation
#               0.54 in E1-E3 and -0.18 in E4), which violates the Kronecker assumption of model (8).
#  Environment 4 is predicted: (a) TLUE, all its records masked; (b) TLUE-TA, only its target masked (the
#  auxiliary trait is observed in E4). Models: UNI, MT (unstructured Sigma_A and Sigma_T), MT-GxEdiag (diagonal
#  Sigma_T), all with Omega = I. Reported: within-environment PC with the phenotype and with the true genetic
#  value, bias of the predicted environment mean and mean absolute error.
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R"); source("R/v4_helpers.R")
args <- commandArgs(TRUE); r1 <- as.integer(args[1]); r2 <- as.integer(args[2])
d0 <- readRDS("data/prepared/EYT_4.rds"); J <- 300; I <- 4
lines <- rownames(d0$G)[1:J]; G <- d0$G[lines, lines]; eg <- eigen(G, symmetric = TRUE); LG <- t(chol(G))
line <- rep(lines, I); env <- rep(paste0("E", 1:I), each = J)
Vg <- c(0.4, 0.7); Rv <- c(0.6, 0.3); Rt <- diag(Rv); Rt[1, 2] <- Rt[2, 1] <- 0.3 * sqrt(prod(Rv))
cc <- function(r) matrix(c(1, r, r, 1), 2, 2)
mc <- list(nIter = 2000, burnIn = 500, thin = 3)
RES <- NULL
for (rep in r1:r2) for (scn in c("constant", "changing")) {
  set.seed(40000 + rep + 1000 * (scn == "changing"))
  Z <- LG %*% matrix(rnorm(J * 2), J, 2) %*% chol(cc(0.3))
  rd <- if (scn == "changing") c(0.9, 0.9, 0.9, -0.9) else rep(0.9, 4)
  g <- matrix(0, J * I, 2)
  for (i in 1:I) { Dd <- LG %*% matrix(rnorm(J * 2), J, 2) %*% chol(cc(rd[i])); g[(i - 1) * J + 1:J, ] <- sqrt(0.6) * Z + sqrt(0.4) * Dd }
  g <- g %*% diag(sqrt(Vg))
  mu <- cbind(rnorm(I, 0, 0.7), rnorm(I, 0, 0.7))
  L <- g + mu[rep(1:I, each = J), ] + matrix(rnorm(J * I * 2), J * I, 2) %*% chol(Rt)
  Y <- data.frame(tgt = 5 + L[, 1], aux = 20 + 2 * L[, 2])
  te <- which(env == "E4")
  ev <- function(p, model, mask) data.frame(rep = rep, scenario = scn, mask = mask, model = model,
                                            PC_pheno = cor(p, Y$tgt[te]), PC_g = cor(p, g[te, 1]), bias = mean(p) - mean(Y$tgt[te]),
                                            MAE = mean(abs(p - Y$tgt[te])))
  Yu <- Y; Yu$tgt[te] <- NA
  set.seed(rep + 5)
  fu <- BMTMEMix(Yu[, "tgt", drop = FALSE], "continuous", line, env, G, covE = "ID", nIter = mc$nIter, burnIn = mc$burnIn, thin = mc$thin, eigG = eg)
  RES <- rbind(RES, ev(fu$pred[fu$cell, 1][te], "UNI", "TLUE"), ev(fu$pred[fu$cell, 1][te], "UNI", "TLUE-TA"))
  for (mask in c("TLUE", "TLUE-TA")) {
    Ym <- Yu; if (mask == "TLUE") Ym$aux[te] <- NA
    for (gx in c("UN", "DIAG")) {
      set.seed(rep + 11 + (gx == "DIAG"))
      fit <- BMTMEMix(Ym, c("continuous", "continuous"), line, env, G, covT = "UN", covTgxe = gx, covE = "ID",
                      nIter = mc$nIter, burnIn = mc$burnIn, thin = mc$thin, eigG = eg)
      RES <- rbind(RES, ev(fit$pred_sys[fit$cell, 1][te], if (gx == "UN") "MT" else "MT-GxEdiag", mask))
    }
  }
  cat("rep", rep, scn, "done\n")
}
dir.create(file.path(V4DIR, "sim_envcor"), showWarnings = FALSE, recursive = TRUE)
fwrite(RES, file.path(V4DIR, "sim_envcor", sprintf("envcor_%03d_%03d.csv", r1, r2)))
cat("done\n")
