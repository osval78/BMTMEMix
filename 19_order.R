# =============================================================================
# 19_order.R — sensitivity of the order-dependent structures to arbitrary and justified orders (v4)
# Usage: Rscript R/19_order.R <dataset> <traits|environments>
#  Scenario 2, TLTE fold 1 (as R/08_sensitivity.R). Structures AR1, ARH1, TOEP, ANTE1, SP.
#  traits      : all 12 distinct orders of the four traits (an order and its reverse give the same set of
#                neighbours and are counted once), in model (8) with unstructured Omega. The "justified" order is
#                the seriation that maximizes the sum of absolute genetic correlations (Sigma_G of the full-data
#                fit) between neighbours.
#  environments: separable model; original order (wheat: file order; rice: chronological), its reverse, four random
#                permutations whose neighbour sets differ from the original one, and a justified order (ascending
#                mean grain yield of the environments, i.e. a stress gradient).
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R"); source("R/v4_helpers.R")
args <- commandArgs(TRUE); d <- args[1]; what <- args[2]
D <- readRDS(sprintf("data/prepared/%s.rds", d)); eg <- eigen(D$G, symmetric = TRUE); types <- D$types_S2
fo <- make_folds("TLTE", D$line, D$env, NFOLD, FOLD_SEED)[[1]]
te <- fo$test; train <- setdiff(seq_len(nrow(D$raw)), te)
Yd <- make_responses(D, train); Ytr <- Yd; Ytr[te, ] <- NA
perms <- function(n) { if (n == 1) return(matrix(1, 1, 1)); p <- perms(n - 1); do.call(rbind, lapply(1:n, function(i) cbind(i, matrix(setdiff(1:n, i)[p], ncol = n - 1)))) }
nbrs <- function(o) paste(sort(apply(cbind(head(o, -1), o[-1]), 1, function(x) paste(sort(x), collapse = "-"))), collapse = "|")
evalfit <- function(fit, ord = seq_along(types)) {
  out <- NULL
  for (k in seq_along(ord)) {
    t <- ord[k]; pr <- fit$pred[fit$cell, k]; pb <- if (types[t] == "ordinal") fit$probs[[k]][fit$cell, , drop = FALSE] else NULL
    for (e in unique(D$env[te])) { idx <- te[D$env[te] == e]
      mt <- metric_trait(types[t], Yd[idx, t], pr[idx], if (is.null(pb)) NULL else pb[idx, , drop = FALSE])
      out <- rbind(out, data.frame(trait = names(Yd)[t], env = e, M1 = mt[["M1"]])) }
  }
  aggregate(M1 ~ trait, out, mean)
}
res <- NULL
if (what == "traits") {
  P <- perms(4); P <- P[!duplicated(apply(P, 1, nbrs)), ]
  SG <- abs(cov2cor(as.matrix(read.csv(sprintf("outputs/tables/SigmaG_%s_S2.csv", d), row.names = 1))))
  score <- apply(P, 1, function(o) sum(SG[cbind(head(o, -1), o[-1])]))
  just <- which.max(score)
  for (st in c("AR1", "ARH1", "TOEP", "ANTE1", "SP")) for (r in seq_len(nrow(P))) {
    ord <- P[r, ]
    set.seed(FIT_SEED + 1000)
    fit <- BMTMEMix(Ytr[, ord], types[ord], D$line, D$env, D$G, covT = st, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                    thin = CV_MCMC$thin, eigG = eg, envLevels = D$envOrder)
    ev <- evalfit(fit, ord)
    res <- rbind(res, data.frame(dataset = d, analysis = "order-traits", structure = st, order = paste(names(Yd)[ord], collapse = "-"),
                                 neighbours = nbrs(ord), justified = r == just, original = all(ord == 1:4) || all(ord == 4:1), ev))
    cat(st, r, "\n")
  }
} else {
  I <- length(D$envOrder); set.seed(99)
  orig <- seq_len(I); cand <- list(original = orig, reversed = rev(orig))
  mgy <- tapply(D$raw[[D$spec[["cont"]]]], D$env, mean)[D$envOrder]; cand$justified <- order(mgy)
  k <- 0
  while (k < 4) { o <- sample(I); if (!nbrs(o) %in% sapply(cand, nbrs)) { k <- k + 1; cand[[paste0("random", k)]] <- o } }
  for (st in c("AR1", "ARH1", "TOEP", "ANTE1", "SP")) for (nm in names(cand)) {
    el <- D$envOrder[cand[[nm]]]
    set.seed(FIT_SEED + 1000)
    fit <- BMTMEMix(Ytr, types, D$line, D$env, D$G, gxe = "separable", covE = st, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                    thin = CV_MCMC$thin, eigG = eg, envLevels = el)
    ev <- evalfit(fit)
    res <- rbind(res, data.frame(dataset = d, analysis = "order-environments", structure = st, order = paste(el, collapse = "-"),
                                 neighbours = nbrs(cand[[nm]]), justified = nm == "justified", original = nm %in% c("original", "reversed"), label = nm, ev))
    cat(st, nm, "\n")
  }
}
dir.create(file.path(V4DIR, "order"), showWarnings = FALSE, recursive = TRUE)
fwrite(res, file.path(V4DIR, "order", sprintf("%s_%s.csv", d, what)))
cat("done", d, what, "\n")
