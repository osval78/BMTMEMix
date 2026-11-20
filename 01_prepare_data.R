# =============================================================================
# 01_prepare_data.R — genomic relationship matrices and mixed-type responses
# for the six data sets (EYT_4, EYT_5, EYT_6, EYT_7, Indica, Japonica).
#   continuous : GY (grain yield)
#   count      : EYT -> DTHD rounded to integer days; rice -> PH rounded to integer cm
#   categorical: EYT -> DTMT; rice -> PHR, in three classes defined by the
#                within-environment tertiles (1 = lower, 2 = middle, 3 = upper third).
#                Scenario 1 treats it as nominal, scenario 2 as ordinal (1 < 2 < 3).
#   binary     : EYT -> Height; rice -> GC, 1 if above the within-environment median.
# =============================================================================
suppressPackageStartupMessages(library(data.table))
source("R/BMTMEMix.R")
dir.create("data/prepared", showWarnings = FALSE)
tertile_class <- function(y, env) {             # within-environment tertiles -> 1, 2, 3
  out <- integer(length(y))
  for (e in unique(env)) {
    i <- env == e; br <- quantile(y[i], c(0, 1/3, 2/3, 1), na.rm = TRUE)
    out[i] <- as.integer(cut(y[i], br, include.lowest = TRUE, labels = FALSE))
  }
  out
}
above_median <- function(y, env) {
  out <- integer(length(y))
  for (e in unique(env)) { i <- env == e; out[i] <- as.integer(y[i] > median(y[i])) }
  out
}
spec <- list(
  EYT_4 = c(cont = "GY", count = "DTHD", cat = "DTMT", bin = "Height"),
  EYT_5 = c(cont = "GY", count = "DTHD", cat = "DTMT", bin = "Height"),
  EYT_6 = c(cont = "GY", count = "DTHD", cat = "DTMT", bin = "Height"),
  EYT_7 = c(cont = "GY", count = "DTHD", cat = "DTMT", bin = "Height"),
  Indica   = c(cont = "GY", count = "PH", cat = "PHR", bin = "GC"),
  Japonica = c(cont = "GY", count = "PH", cat = "PHR", bin = "GC"))
summ <- NULL
for (d in names(spec)) {
  X <- fread(sprintf("data/X_%s.csv", d)); Y <- as.data.frame(fread(sprintf("data/Y_%s.csv", d)))
  if (names(X)[1] %in% c("V1", "")) X <- X[, -1]
  if (names(Y)[1] %in% c("V1", "")) Y <- Y[, -1]
  ids <- X$Line; M <- as.matrix(X[, -1]); rownames(M) <- ids
  G <- make_G(M); rownames(G) <- colnames(G) <- ids
  env <- as.character(Y$Env); line <- as.character(Y$Line)
  s <- spec[[d]]
  Yd <- data.frame(Y[[s["cont"]]], round(Y[[s["count"]]]), tertile_class(Y[[s["cat"]]], env), above_median(Y[[s["bin"]]], env))
  names(Yd) <- unname(s)
  nu <- attr(G, "n_used"); G <- G[unique(line), unique(line)]
  saveRDS(list(Yd = Yd, types = c("continuous", "count", "nominal", "binary"),
               types_S2 = c("continuous", "count", "ordinal", "binary"),
               line = line, env = env, G = G, raw = Y, spec = s, envOrder = sort(unique(env)),
               n_markers = ncol(M), n_markers_used = nu),
          sprintf("data/prepared/%s.rds", d))
  tab <- table(Yd[[3]]); bt <- table(Yd[[4]])
  summ <- rbind(summ, data.frame(dataset = d, J = nrow(G), I = length(unique(env)), n = nrow(Yd),
                                 cells = nrow(G) * length(unique(env)), markers = ncol(M), markers_used = nu,
                                 envs = paste(sort(unique(env)), collapse = ", "),
                                 continuous = s["cont"], count = s["count"], categorical = s["cat"], binary = s["bin"],
                                 count_range = paste(range(Yd[[2]]), collapse = "-"),
                                 cat_freq = paste(tab, collapse = "/"), bin_freq = paste(bt, collapse = "/")))
  cat(d, ": J =", nrow(G), " I =", length(unique(env)), " n =", nrow(Yd), " markers =", ncol(M), "\n")
}
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(summ, "outputs/tables/datasets_summary.csv", row.names = FALSE)
print(summ)
