# =============================================================================
# 06_summarize.R — summary tables (outputs/tables) and figures (figures/) — version 4
# v4: every accuracy reported in the tables is the pooled-within-environment estimand (metric computed within each
#     environment on the pooled out-of-fold predictions, averaged over environments), the same one used for the
#     paired bootstrap of the differences; 1,000 bootstrap replicates (BMTME_BOOT); bootstraps run in parallel.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(parallel) })
NC <- as.integer(Sys.getenv("BMTME_CORES", "2"))
source("R/00_settings.R")
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE); dir.create("figures", showWarnings = FALSE)
set.seed(2026)
COL <- c(BMTMEMix = "#2a78d6", Univariate = "#eb6834")
theme_paper <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
        strip.text = element_text(face = "bold"), legend.position = "top",
        plot.background = element_rect(fill = "white", colour = NA))
DS <- DATASETS; SCH <- c(SCHEMES, "CHRONO")
ROLE_LAB <- c(continuous = "Continuous", count = "Count", categorical = "Categorical", binary = "Binary")
role_of <- function(type) ifelse(type %in% c("nominal", "ordinal"), "categorical", type)
rdcsv <- function(pattern, dir) { f <- list.files(dir, pattern = pattern, full.names = TRUE); if (!length(f)) return(NULL); rbindlist(lapply(f, fread), fill = TRUE) }
rdpreds <- function(dir) { f <- list.files(dir, pattern = "\\.rds$", full.names = TRUE); if (!length(f)) return(NULL)
  rbindlist(lapply(f, function(x) as.data.table(readRDS(x)$preds)), fill = TRUE) }
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))

## =============================== 1. cross-validation: fold means ===============================
S1 <- rdcsv("_metrics.csv$", "outputs/cv/S1"); S2 <- rdcsv("_metrics.csv$", "outputs/cv/S2")
reuse <- S1[model == "Univariate" & type %in% c("continuous", "count", "binary")][, scenario := "S2"]
M <- rbind(S1, S2, reuse, fill = TRUE)[, role := role_of(type)]
fwrite(M, "outputs/tables/cv_metrics_all_folds_envs.csv")
F1 <- M[, .(M1 = mean(M1, na.rm = TRUE), M2 = mean(M2, na.rm = TRUE), M3 = mean(M3, na.rm = TRUE), MAE = mean(MAE, na.rm = TRUE), secs = mean(secs)),
        by = .(dataset, scenario, scheme, trait, type, role, model, fold)]
S <- F1[, .(M1 = mean(M1), M1_se = se(M1), M2 = mean(M2), M2_se = se(M2), M3 = mean(M3), MAE = mean(MAE), nfold = .N, secs = mean(secs)),
        by = .(dataset, scenario, scheme, trait, type, role, model)]
S[, metric1 := ifelse(role %in% c("continuous", "count"), "PC", "PCCC")][, metric2 := ifelse(role %in% c("continuous", "count"), "NRMSE", "Brier")]
fwrite(S, "outputs/tables/cv_summary.csv")
fwrite(S[, .(secs = mean(secs)), by = .(dataset, scenario, model)], "outputs/tables/cv_seconds_per_fit.csv")

## =============================== 2. record-level predictions ===============================
P1 <- rdpreds("outputs/cv/S1"); P2 <- rdpreds("outputs/cv/S2")
P2 <- rbind(P2, P1[model == "Univariate" & type %in% c("continuous", "count", "binary")][, scenario := "S2"], fill = TRUE)
P <- rbind(P1, P2, fill = TRUE)[, role := role_of(type)]
fwrite(P, "outputs/tables/cv_predictions_all.csv.gz")   # record-level predictions (compressed)
# metrics of models A and B on a set of records, computed within environments and averaged over environments
mvals <- function(o, p, pr, role) {
  if (role %in% c("continuous", "count")) c(suppressWarnings(cor(o, p)), sqrt(mean((o - p)^2)) / mean(o), NA_real_)
  else if (role == "binary") c(mean((p > 0.5) == (o == 1)), mean((p - o)^2), NA_real_)
  else { Yi <- cbind(o == 1, o == 2, o == 3) * 1
    c(mean(max.col(pr, ties.method = "first") == o), mean(rowSums((pr - Yi)^2)),
      mean(((pr[, 1] - Yi[, 1])^2 + (pr[, 1] + pr[, 2] - Yi[, 1] - Yi[, 2])^2) / 2)) }
}
envmetrics <- function(w, rows = NULL) {
  if (!is.null(rows)) w <- w[rows]
  rl <- w$role[1]
  v <- sapply(split(seq_len(nrow(w)), w$env), function(ix) {
    x <- w[ix]
    c(mvals(x$obs, x$A, cbind(x$A_p1, x$A_p2, x$A_p3), rl), mvals(x$obs, x$B, cbind(x$B_p1, x$B_p2, x$B_p3), rl)) })
  rowMeans(v, na.rm = TRUE)
}
# wide table: one row per record, columns for model A and model B
widen <- function(dt, A, B) {
  a <- dt[model == A, .(fold, line, env, obs, role, A = pred, A_p1 = p1, A_p2 = p2, A_p3 = p3)]
  b <- dt[model == B, .(fold, line, env, B = pred, B_p1 = p1, B_p2 = p2, B_p3 = p3)]
  merge(a, b, by = c("fold", "line", "env"))
}
boot_diff <- function(w, B = as.integer(Sys.getenv("BMTME_BOOT", "1000"))) {
  # paired cluster bootstrap over lines (all records of a resampled line enter together, in every environment)
  w <- as.data.table(w); lines <- unique(w$line); idx <- split(seq_len(nrow(w)), w$line)
  e0 <- envmetrics(w); est <- e0[4:6] - e0[1:3]
  bs <- replicate(B, { e <- envmetrics(w, unlist(idx[sample(lines, replace = TRUE)], use.names = FALSE)); e[4:6] - e[1:3] })
  data.table(metric = c("M1", "M2", "M3"), A = e0[1:3], B = e0[4:6], diff = est,
             lo = apply(bs, 1, quantile, 0.025, na.rm = TRUE), hi = apply(bs, 1, quantile, 0.975, na.rm = TRUE))
}
jobs <- list()
for (sc in c("S1", "S2")) for (d in DS) for (s in SCH) for (tr in unique(P[dataset == d, trait])) {
  dt <- P[scenario == sc & dataset == d & scheme == s & trait == tr]
  if (!nrow(dt) || !all(c("BMTMEMix", "Univariate") %in% dt$model)) next
  jobs[[length(jobs) + 1]] <- list(w = widen(dt, "Univariate", "BMTMEMix"), sc = sc, d = d, s = s, tr = tr)
}
BD <- rbindlist(mclapply(seq_along(jobs), function(k) { j <- jobs[[k]]; set.seed(2026 + k)
  b <- boot_diff(j$w); b[, `:=`(scenario = j$sc, dataset = j$d, scheme = j$s, trait = j$tr, role = j$w$role[1])]; b }, mc.cores = NC))
BD <- BD[!(role %in% c("continuous", "count", "binary") & metric == "M3")]
fwrite(BD, "outputs/tables/diff_MT_minus_UNI_bootstrap.csv")
# pooled-within-environment accuracies (A = univariate, B = BMTMEMix): the estimand of every table in v4
fwrite(BD[, .(scenario, dataset, scheme, trait, role, metric, UNI = A, MT = B, diff, lo, hi)], "outputs/tables/cv_pooled_summary.csv")
# descriptive summary across data sets (heterogeneity shown by the SD; CIs counted per data set)
BDs <- BD[, .(mean_diff = mean(diff), sd_diff = sd(diff), n = .N,
              n_better = sum((metric == "M1" & lo > 0) | (metric != "M1" & hi < 0)),
              n_worse = sum((metric == "M1" & hi < 0) | (metric != "M1" & lo > 0))), by = .(scenario, scheme, role, metric)]
fwrite(BDs, "outputs/tables/diff_MT_minus_UNI_summary.csv")
fwrite(BD[, .(mean_diff = mean(diff), sd_diff = sd(diff), n = .N,
              n_better = sum((metric == "M1" & lo > 0) | (metric != "M1" & hi < 0)),
              n_worse = sum((metric == "M1" & hi < 0) | (metric != "M1" & lo > 0))), by = .(scenario, metric)], "outputs/tables/diff_MT_minus_UNI_overall.csv")

## =============================== 3. nominal vs ordinal ===============================
jobs <- list()
for (m in c("BMTMEMix", "Univariate")) for (d in DS) for (s in SCH) {
  a <- P[scenario == "S1" & model == m & dataset == d & scheme == s & role == "categorical"]
  b <- P[scenario == "S2" & model == m & dataset == d & scheme == s & role == "categorical"]
  if (!nrow(a) || !nrow(b)) next
  dt <- rbind(a[, model := "nominal"], b[, model := "ordinal"])
  jobs[[length(jobs) + 1]] <- list(w = widen(dt, "nominal", "ordinal"), m = m, d = d, s = s)
}
SC <- rbindlist(mclapply(seq_along(jobs), function(k) { j <- jobs[[k]]; set.seed(3026 + k)
  bb <- boot_diff(j$w); bb[, `:=`(model = j$m, dataset = j$d, scheme = j$s)]; bb }, mc.cores = NC))
if (!nrow(SC)) SC <- NULL
if (is.null(SC)) SC <- data.table(metric = character(), A = numeric(), B = numeric(), diff = numeric(), lo = numeric(), hi = numeric(), model = character(), dataset = character(), scheme = character())
fwrite(SC, "outputs/tables/diff_ordinal_minus_nominal_bootstrap.csv")
fwrite(SC[, .(mean_nominal = mean(A), mean_ordinal = mean(B), mean_diff = mean(diff), sd_diff = sd(diff), n = .N,
              n_ord_better = sum((metric == "M1" & lo > 0) | (metric != "M1" & hi < 0)),
              n_ord_worse = sum((metric == "M1" & hi < 0) | (metric != "M1" & lo > 0))), by = .(model, scheme, metric)],
       "outputs/tables/diff_ordinal_minus_nominal_summary.csv")
fwrite(SC[, .(mean_nominal = mean(A), mean_ordinal = mean(B), mean_diff = mean(diff), n = .N,
              n_ord_better = sum((metric == "M1" & lo > 0) | (metric != "M1" & hi < 0)),
              n_ord_worse = sum((metric == "M1" & hi < 0) | (metric != "M1" & lo > 0))), by = .(model, metric)],
       "outputs/tables/diff_ordinal_minus_nominal_overall.csv")
# effect of the scenario on the other traits (BMTMEMix)
O <- S[model == "BMTMEMix" & role != "categorical"]
OW <- dcast(O, dataset + scheme + trait + role ~ scenario, value.var = "M1")
if ("S2" %in% names(OW)) fwrite(OW[, .(mean_abs_change = mean(abs(S2 - S1)), max_abs_change = max(abs(S2 - S1)), mean_change = mean(S2 - S1)), by = role],
       "outputs/tables/scenario_other_traits_by_role.csv")

## =============================== 4. calibration of probabilistic predictions ===============================
events <- function(dt) {          # binary event and cumulative ordinal/nominal events: (event, p, y)
  rbind(dt[role == "binary", .(event = "Binary: P(Y = 1)", p = pred, y = obs)],
        dt[role == "categorical", .(event = "Categorical: P(Y <= 1)", p = p1, y = as.numeric(obs <= 1))],
        dt[role == "categorical", .(event = "Categorical: P(Y <= 2)", p = p1 + p2, y = as.numeric(obs <= 2))])
}
CAL <- NULL; REL <- NULL
for (sc in c("S1", "S2")) for (m in c("BMTMEMix", "Univariate")) for (s in SCH) {
  dt <- P[scenario == sc & model == m & scheme == s]; if (!nrow(dt)) next
  ev <- events(dt)
  for (e in unique(ev$event)) {
    x <- ev[event == e & is.finite(p)]; x[, pp := pmin(pmax(p, 1e-4), 1 - 1e-4)]
    fitc <- tryCatch(coef(glm(y ~ qlogis(pp), family = binomial, data = x)), error = function(err) c(NA, NA))
    x[, bin := cut(pp, breaks = seq(0, 1, 0.1), include.lowest = TRUE)]
    rb <- x[, .(mean_p = mean(pp), obs_freq = mean(y), n = .N), by = bin]
    CAL <- rbind(CAL, data.table(scenario = sc, model = m, scheme = s, event = e, cal_intercept = fitc[1], cal_slope = fitc[2],
                                 ece = sum(rb$n * abs(rb$mean_p - rb$obs_freq)) / sum(rb$n), n = nrow(x)))
    REL <- rbind(REL, rb[, `:=`(scenario = sc, model = m, scheme = s, event = e)])
  }
}
fwrite(CAL, "outputs/tables/calibration_summary.csv"); fwrite(REL, "outputs/tables/reliability_bins.csv")
REL[, grp := paste(model, ifelse(scenario == "S1", "nominal", "ordinal"))]
pr <- REL[scheme %in% c("TLTE", "TLUE") & n >= 20]
p <- ggplot(pr, aes(mean_p, obs_freq, colour = grp)) + geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.2) + facet_grid(scheme ~ event) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = c("BMTMEMix nominal" = "#86b6ef", "BMTMEMix ordinal" = "#1c5cab", "Univariate nominal" = "#f5b08f", "Univariate ordinal" = "#c24a1a"), name = NULL) +
  labs(x = "Mean predicted probability (decile bins)", y = "Observed frequency") + theme_paper
ggsave("figures/Fig_reliability.png", p, width = 7.4, height = 5.2, dpi = 300)

## =============================== 5. untested environments: absolute errors, calibration, intervals ===============================
UE <- P[scheme %in% c("TLUE", "ULUE", "CHRONO") & role %in% c("continuous", "count")]
UEs <- UE[, .(PC = suppressWarnings(cor(obs, pred)), bias = mean(pred - obs), rel_bias = mean(pred - obs) / mean(obs),
              MAE = mean(abs(pred - obs)), slope = if (sd(pred) > 0) coef(lm(obs ~ pred))[2] else NA_real_,
              cover90 = mean(abs(obs - pred) <= 1.645 * pred_sd, na.rm = TRUE), n = .N),
          by = .(scenario, dataset, scheme, env, trait, role, model)]
fwrite(UEs, "outputs/tables/untested_env_by_environment.csv")
fwrite(UEs[, .(PC = mean(PC), abs_rel_bias = mean(abs(rel_bias)), slope = median(slope, na.rm = TRUE), cover90 = mean(cover90, na.rm = TRUE), n_env = .N),
           by = .(scenario, scheme, role, model)], "outputs/tables/untested_env_summary.csv")

# posterior predictive check of the count model (TLTE): observed vs predicted within-environment dispersion
# (predicted variance = variance of the predictive means + mean predictive variance), share of zeros, 90% PI coverage
PPC <- P[scheme == "TLTE" & role == "count" & model == "BMTMEMix" & scenario == "S1",
         .(var_obs = var(obs), var_pred = var(pred) + mean(pred_sd^2), mean_obs = mean(obs), mean_pred = mean(pred),
           zeros_obs = mean(obs == 0), cover90 = mean(abs(obs - pred) <= 1.645 * pred_sd), n = .N), by = .(dataset, env)]
PPC[, disp_ratio := var_pred / var_obs]
fwrite(PPC, "outputs/tables/ppc_count.csv")

## =============================== 6. trait-assisted, competitors ===============================
TA <- rdcsv("_ULTE-TA_metrics.csv$", "outputs/cv_TA"); TAX <- rdcsv("_ULTE-TA-xenv_metrics.csv$", "outputs/cv_TA")
if (!is.null(TA)) {
TAC <- rdcsv("_ULTE-TA-competitor_metrics.csv$", "outputs/cv_TA")   # nested-CV competitor (R/07b_ta_competitor.R)
TA <- TA[model != "UNI+secondaries"]                                  # in-sample version written by 07 is not used
TT <- rbind(TA, TAX, TAC, fill = TRUE)[, role := role_of(type)]
T1 <- TT[, .(M1 = mean(M1, na.rm = TRUE), M2 = mean(M2, na.rm = TRUE), M3 = mean(M3, na.rm = TRUE)), by = .(dataset, scenario, scheme, trait, role, model, fold)]
Ts <- T1[, .(M1 = mean(M1), M1_se = se(M1), M2 = mean(M2), M3 = mean(M3)), by = .(dataset, scenario, scheme, trait, role, model)]
ul <- S[scheme == "ULTE", .(dataset, scenario, trait, role, model = paste0(model, " (ULTE, all traits masked)"), M1, M1_se, M2, M3, scheme = "ULTE")]
TAall <- rbind(Ts, ul, fill = TRUE)
fwrite(TAall, "outputs/tables/trait_assisted_all.csv")
TW <- dcast(TAall[scheme %in% c("ULTE", "ULTE-TA")], dataset + scenario + trait + role ~ model, value.var = "M1")
fwrite(TW, "outputs/tables/trait_assisted_wide.csv")
if (any(TAall$scheme == "ULTE-TA-xenv")) TX <- dcast(TAall[scheme %in% c("ULTE-TA-xenv")], dataset + scenario + trait + role ~ model, value.var = "M1")
if (exists("TX")) fwrite(TX, "outputs/tables/trait_assisted_xenv_wide.csv")
p <- ggplot(TAall[scheme %in% c("ULTE", "ULTE-TA") & model %in% c("Univariate (ULTE, all traits masked)", "BMTMEMix (ULTE, all traits masked)",
                                                         "UNI+secondaries", "BMTMEMix (phenotypic)", "BMTMEMix (systematic)")],
            aes(scenario, M1, fill = factor(model, c("Univariate (ULTE, all traits masked)", "BMTMEMix (ULTE, all traits masked)", "UNI+secondaries",
                                                    "BMTMEMix (systematic)", "BMTMEMix (phenotypic)")))) +
  geom_col(position = position_dodge(width = 0.86), width = 0.82, colour = "white", linewidth = 0.2) +
  facet_grid(factor(ROLE_LAB[role], ROLE_LAB) ~ factor(dataset, DS), scales = "free_y") +
  scale_fill_manual(values = c("#eb6834", "#86b6ef", "#eda100", "#5598e7", "#1c5cab"), name = NULL,
                    labels = c("UNI", "MT, all traits masked", "UNI + secondaries (regression)", "MT-TA systematic", "MT-TA phenotypic")) +
  labs(x = "Response scenario", y = "Accuracy (PC or PCCC), untested lines") + theme_paper + theme(legend.text = element_text(size = 7))
ggsave("figures/Fig_trait_assisted.png", p, width = 7.6, height = 6.6, dpi = 300)
}

## =============================== 7. other competitors and variants ===============================
G <- rdcsv("_metrics.csv$", "outputs/cv/S1_gauss")
if (!is.null(G)) {
  g1 <- G[, .(M1 = mean(M1), M2 = mean(M2), M3 = mean(M3, na.rm = TRUE)), by = .(dataset, scheme, trait, type, fold)][, .(M1_gauss = mean(M1), M2_gauss = mean(M2), M3_gauss = mean(M3)), by = .(dataset, scheme, trait)]
  b1 <- S[model == "BMTMEMix" & scheme %in% c("TLTE", "ULTE"), .(dataset, scheme, scenario, trait, role, M1, M2, M3)]
  GW <- merge(b1, g1, by = c("dataset", "scheme", "trait"))
  fwrite(GW, "outputs/tables/gaussian_vs_BMTMEMix.csv")
  fwrite(GW[, .(M1_mix = mean(M1), M1_gauss = mean(M1_gauss), M2_mix = mean(M2), M2_gauss = mean(M2_gauss), M3_mix = mean(M3), M3_gauss = mean(M3_gauss),
                mix_better_M2 = sum(M2 < M2_gauss), n = .N), by = .(scenario, role)], "outputs/tables/gaussian_vs_BMTMEMix_summary.csv")
}
MG <- rdcsv("_metrics.csv$", "outputs/cv/S2_mcmcglmm")
if (!is.null(MG)) {
  m1 <- MG[, .(M1 = mean(M1), M2 = mean(M2), M3 = mean(M3, na.rm = TRUE), secs = mean(secs)), by = .(dataset, trait, type, fold)]
  b1 <- M[scenario == "S2" & scheme == "TLTE" & dataset %in% unique(MG$dataset) & fold %in% unique(MG$fold)]
  b1 <- b1[, .(M1 = mean(M1), M2 = mean(M2), M3 = mean(M3, na.rm = TRUE), secs = mean(secs)), by = .(dataset, trait, type, model, fold)]
  MGt <- rbind(m1[, model := "MCMCglmm"], b1, fill = TRUE)[, .(M1 = mean(M1), M2 = mean(M2), M3 = mean(M3), secs = mean(secs)), by = .(dataset, trait, type, model)]
  fwrite(MGt, "outputs/tables/mcmcglmm_comparison.csv")
}
SP <- rdcsv("_metrics.csv$", "outputs/cv/S1_sep")
if (!is.null(SP)) {
  a <- SP[, .(M1 = mean(M1)), by = .(dataset, trait, type, fold)][, .(M1_sep = mean(M1)), by = .(dataset, trait, type)]
  b <- S[scenario == "S1" & scheme == "TLTE" & model == "BMTMEMix", .(dataset, trait, M1_main = M1)]
  fwrite(merge(a, b, by = c("dataset", "trait"))[, dM1 := M1_main - M1_sep], "outputs/tables/separable_vs_main.csv")
}
RF <- rdcsv("_metrics.csv$", "outputs/cv/S1_ref")
if (!is.null(RF)) {
  a <- RF[, .(M1 = mean(M1), M2 = mean(M2), M3 = mean(M3, na.rm = TRUE)), by = .(dataset, trait, type, fold)][, .(M1_ref1 = mean(M1), M2_ref1 = mean(M2), M3_ref1 = mean(M3)), by = .(dataset, trait, type)]
  b <- S[scenario == "S1" & scheme == "TLTE" & model == "BMTMEMix", .(dataset, trait, M1_ref3 = M1, M2_ref3 = M2, M3_ref3 = M3)]
  fwrite(merge(a, b, by = c("dataset", "trait")), "outputs/tables/nominal_reference_sensitivity.csv")
}
UN <- rdcsv("_metrics.csv$", "outputs/cv/S1_unE")
if (!is.null(UN)) {
  a <- UN[, .(M1 = mean(M1)), by = .(dataset, scheme, trait, model, fold)][, .(M1_UN = mean(M1)), by = .(dataset, scheme, trait, model)]
  b <- S[scenario == "S1", .(dataset, scheme, trait, model, M1_ID = M1)]
  UNm <- merge(a, b, by = c("dataset", "scheme", "trait", "model"))
  fwrite(UNm, "outputs/tables/untested_env_UN_vs_ID.csv")
  fwrite(UNm[, .(M1_UN = mean(M1_UN), M1_ID = mean(M1_ID)), by = .(scheme, model)], "outputs/tables/untested_env_UN_vs_ID_summary.csv")
}

## =============================== 8. covariance structures, sensitivity, long chains ===============================
ST <- rdcsv("_structures.csv$", "outputs/structures")
if (!is.null(ST)) {
  ST[, structure := ifelse(target == "traits", covT, covE)]
  A <- ST[, .(M1 = mean(M1), M2 = mean(M2), secs = mean(secs)), by = .(dataset, fold, target, structure, trait, type)]
  fwrite(A, "outputs/tables/structures_by_trait.csv")
  A2 <- A[, .(M1 = mean(M1), secs = mean(secs)), by = .(dataset, fold, target, structure)][, .(M1 = mean(M1), M1_fold_sd = sd(M1), secs = mean(secs)), by = .(dataset, target, structure)]
  ref <- A2[structure == "UN", .(dataset, target, M1_UN = M1)]
  A2 <- merge(A2, ref, by = c("dataset", "target"))[, dUN := 100 * (M1 - M1_UN) / M1_UN]
  fwrite(A2, "outputs/tables/structures_mean.csv")
  lev <- names(c(UN = 1, DIAG = 2, FA = 3, REC = 4, CS = 5, ID = 6, CSH = 7, AR1 = 8, ARH1 = 9, TOEP = 10, ANTE1 = 11, SP = 12))
  A2[, structure := factor(structure, rev(lev))][, dataset := factor(dataset, DS)]
  A2[, target := factor(target, c("environments", "traits"), c("Structure for environments (separable model)", "Structure for traits (main + G×E model)"))]
  p <- ggplot(A2, aes(dataset, structure, fill = pmax(pmin(dUN, 5), -5))) + geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", M1)), size = 2.3) + facet_wrap(~target) +
    scale_fill_gradient2(low = "#eb6834", mid = "#f0efec", high = "#2a78d6", midpoint = 0, limits = c(-5, 5), name = "Change vs UN (%)") +
    labs(x = NULL, y = NULL) + theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1), plot.background = element_rect(fill = "white", colour = NA))
  ggsave("figures/Fig_structures.png", p, width = 7.6, height = 4.6, dpi = 300)
}
SE <- rdcsv("_sensitivity.csv$", "outputs/sensitivity")
if (!is.null(SE)) fwrite(SE[, .(M1 = mean(M1), M3 = mean(M3, na.rm = TRUE), corG_12 = mean(corG_12), corG_23 = mean(corG_23), corA_34 = mean(corA_34)),
                            by = .(dataset, analysis, setting, trait)], "outputs/tables/sensitivity_summary.csv")
LC <- rdcsv("_long_chain.csv$", "outputs/long_chain")
if (!is.null(LC)) {
  lc <- LC[, .(M1 = mean(M1), M2 = mean(M2), M3 = mean(M3, na.rm = TRUE), secs = mean(secs)), by = .(dataset, fold, nIter, trait, type)]
  lcw <- dcast(lc, dataset + fold + trait + type ~ nIter, value.var = c("M1", "M2", "secs"))
  fwrite(lcw, "outputs/tables/long_chain_stability.csv")
}

## =============================== 9. full fits: diagnostics and identified covariances ===============================
DG <- rdcsv("_diagnostics.csv$", "outputs/fullfit"); PA <- rdcsv("_predagree.csv$", "outputs/fullfit")
if (!is.null(DG)) {
  DG[, family := sub("\\[.*", "", parameter)]
  fwrite(DG, "outputs/tables/convergence_all.csv")
  fwrite(DG[, .(max_rhat = max(rhat), min_ess_bulk = min(ess_bulk), min_ess_tail = min(ess_tail), n = .N), by = .(dataset, scenario, family)], "outputs/tables/convergence_by_family.csv")
  fwrite(DG[, .(max_rhat = max(rhat), median_rhat = median(rhat), min_ess_bulk = min(ess_bulk), median_ess_bulk = median(ess_bulk), min_ess_tail = min(ess_tail),
                share_rhat_ok = mean(rhat < 1.01), n = .N), by = .(dataset, scenario)], "outputs/tables/convergence_summary.csv")
  if (!is.null(PA)) { PA2 <- merge(PA, DG[, .(min_ess_bulk = min(ess_bulk)), by = .(dataset, scenario)], by = c("dataset", "scenario"))
    PA2[, ess_per_second := min_ess_bulk / (secs_per_chain * 3)]; fwrite(PA2, "outputs/tables/chain_agreement_and_efficiency.csv") }
}
DGL <- rdcsv("_diagnostics.csv$", "outputs/fullfit_long")
if (!is.null(DGL) && !is.null(DG)) {
  DGL[, family := sub("\\[.*", "", parameter)]
  LL <- merge(DG[, .(dataset, scenario, parameter, family, mean_6k = mean, rhat_6k = rhat, ess_bulk_6k = ess_bulk)],
              DGL[, .(dataset, scenario, parameter, mean_30k = mean, rhat_30k = rhat, ess_bulk_30k = ess_bulk, ess_tail_30k = ess_tail)],
              by = c("dataset", "scenario", "parameter"))
  fwrite(LL, "outputs/tables/convergence_long_chains.csv")
}
cl <- function(Cm, nm) { d <- as.data.frame(as.table(Cm)); names(d) <- c("r", "c", "v"); d$panel <- nm; d }
VC <- NULL; THR <- NULL; KT <- NULL
for (d in DS) for (sc in c("S1", "S2")) {
  f <- sprintf("outputs/fullfit/%s_%s_fullfit.rds", d, sc); if (!file.exists(f)) next
  fit <- readRDS(f)
  for (nm in c("SigmaA", "SigmaG", "SigmaT", "R")) fwrite(as.data.table(round(fit[[nm]], 6), keep.rownames = "dim"), sprintf("outputs/tables/%s_%s_%s.csv", nm, d, sc))
  VC <- rbind(VC, data.table(dataset = d, scenario = sc, dim = fit$latnames, SigmaA = diag(fit$SigmaA), SigmaG = diag(fit$SigmaG), SigmaTn = diag(fit$SigmaT), R = diag(fit$R)))
  th <- fit$thresholds[[which(fit$types == "ordinal")[1]]]; if (length(th)) THR <- rbind(THR, data.table(dataset = d, scenario = sc, gamma2 = th[2]))
  # total genetic correlations between environments, per trait (identified)
  st <- cumsum(c(1, fit$nd))[seq_along(fit$types)]
  for (k in seq_along(st)) { Kc <- cov2cor(fit$Ktot[, , st[k]]); KT <- rbind(KT, data.table(dataset = d, scenario = sc, trait = fit$latnames[st[k]],
                                   min_r = min(Kc[upper.tri(Kc)]), mean_r = mean(Kc[upper.tri(Kc)]), max_r = max(Kc[upper.tri(Kc)]))) }
  Kgy <- cov2cor(fit$Ktot[, , 1]); dimnames(Kgy) <- list(fit$envs, fit$envs)
  D <- rbind(cl(cov2cor(fit$SigmaA), "a) Main-effect genetic correlations (Sigma_A)"), cl(cov2cor(fit$SigmaG), "b) Within-environment genetic correlations (Sigma_G)"),
             cl(cov2cor(fit$R), "c) Residual correlations (R)"), cl(Kgy, paste0("d) Total genetic correlations between environments, ", fit$latnames[1])))
  p <- ggplot(D, aes(c, r, fill = v)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.2f", v)), size = 2.4) +
    facet_wrap(~panel, scales = "free", ncol = 2) + scale_fill_gradient2(low = "#eb6834", mid = "#f0efec", high = "#2a78d6", limits = c(-1, 1), name = "Correlation") +
    scale_y_discrete(limits = rev) + labs(x = NULL, y = NULL) + theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold", hjust = 0), plot.background = element_rect(fill = "white", colour = NA))
  ggsave(sprintf("figures/Fig_correlations_%s_%s.png", d, sc), p, width = 7.2, height = 6.2, dpi = 250)
  # trace plots of the three chains for selected identified parameters
  dr <- fit$draws; keep <- intersect(names(dr[[1]]), c(grep("^SigmaG\\[", names(dr[[1]]), value = TRUE), "gamma_2", grep("^R\\[", names(dr[[1]]), value = TRUE)))
  keep <- keep[seq_len(min(9, length(keep)))]
  TL <- rbindlist(lapply(seq_along(dr), function(ch) data.table(chain = factor(ch), iter = seq_len(nrow(dr[[ch]])), stack(dr[[ch]][, keep, drop = FALSE]))))
  p <- ggplot(TL, aes(iter, values, colour = chain)) + geom_line(linewidth = 0.2) + facet_wrap(~ind, scales = "free_y", ncol = 3) +
    scale_colour_manual(values = c("#2a78d6", "#eb6834", "#1baf7a")) + labs(x = "Saved draw", y = NULL, title = paste("Three dispersed chains —", d, sc)) + theme_paper
  ggsave(sprintf("figures/FigS_trace_%s_%s.png", d, sc), p, width = 7.2, height = 5, dpi = 150)
}
if (!is.null(VC)) fwrite(VC, "outputs/tables/varcomp_all.csv")
if (!is.null(THR)) fwrite(THR, "outputs/tables/thresholds_fullfit.csv")
if (!is.null(KT)) fwrite(KT, "outputs/tables/total_genetic_env_correlations.csv")

## =============================== 10. simulation and SBC ===============================
SPm <- rdcsv("^params_", "outputs/simulation"); SAc <- rdcsv("^accuracy_", "outputs/simulation"); GV <- rdcsv("^gv_", "outputs/simulation")
if (!is.null(SPm)) {
  SPm[, `:=`(err = mean - true, cover = (true >= q05 & true <= q95))]
  fwrite(SPm, "outputs/tables/simulation_parameters_all.csv")
  fwrite(SPm[, .(bias = mean(err), rmse = sqrt(mean(err^2)), coverage90 = mean(cover), n_par = uniqueN(parameter), n_rep = uniqueN(rep)), by = .(setting, family)],
         "outputs/tables/simulation_by_family.csv")
  fwrite(SPm[, .(true = mean(true), mean_est = mean(mean), bias = mean(err), rmse = sqrt(mean(err^2)), coverage90 = mean(cover), n_rep = .N), by = .(setting, family, parameter)],
         "outputs/tables/simulation_by_parameter.csv")
  fwrite(unique(SPm[, .(setting, rep, freq, secs)]), "outputs/tables/simulation_replicates.csv")
}
if (!is.null(SAc)) fwrite(SAc[, .(mean_cor = mean(cor_g), sd_cor = sd(cor_g), n = .N), by = .(setting, target, model)], "outputs/tables/simulation_genetic_accuracy.csv")
if (!is.null(GV)) fwrite(GV[, .(mean_cor = mean(cor_g_full), sd = sd(cor_g_full), n = .N), by = .(setting, dim)], "outputs/tables/simulation_gv_full.csv")
SB <- rbind(rdcsv("^sbc_", "outputs/simulation")[, sampler := "basic"], rdcsv("^sbcpx_", "outputs/simulation")[, sampler := "with scale move"], fill = TRUE)
if (!is.null(SB) && nrow(SB)) {
  chi <- SB[, { tab <- table(cut(rank, breaks = seq(-0.5, 99.5, length.out = 11))); ct <- chisq.test(tab); .(chisq = unname(ct$statistic), p_value = ct$p.value, n = .N) }, by = .(sampler, parameter)]
  fwrite(chi, "outputs/tables/sbc_uniformity.csv")
  hl <- SB[, .(h = .N / 10), by = .(sampler, parameter)]
  p <- ggplot(SB, aes(rank)) + geom_histogram(breaks = seq(-0.5, 99.5, length.out = 11), fill = "#2a78d6", colour = "white") +
    geom_hline(data = hl, aes(yintercept = h), colour = "#eb6834", linewidth = 0.4) + facet_grid(sampler ~ parameter, scales = "free_y") +
    labs(x = "Rank of the true value among 99 posterior draws", y = "Replicates") + theme_paper + theme(strip.text = element_text(size = 6), axis.text = element_text(size = 6))
  ggsave("figures/FigS_sbc.png", p, width = 9, height = 3.6, dpi = 200)
}

## =============================== 11. main figures of cross-validation ===============================
PS <- rbind(BD[metric == "M1", .(scenario, dataset, scheme, role, model = "Univariate", M1 = A)], BD[metric == "M1", .(scenario, dataset, scheme, role, model = "BMTMEMix", M1 = B)])
PS[, rlab := factor(ROLE_LAB[role], ROLE_LAB)]
for (sc in c("S1", "S2")) {
  d <- PS[scenario == sc & scheme %in% SCHEMES]; d[, scheme := factor(scheme, SCHEMES)]; d[, dataset := factor(dataset, DS)]
  p <- ggplot(d, aes(scheme, M1, fill = model)) + geom_col(position = position_dodge(width = 0.8), width = 0.72, colour = "white", linewidth = 0.3) +
    facet_grid(rlab ~ dataset, scales = "free_y") + scale_fill_manual(values = COL, name = NULL) +
    labs(x = "Cross-validation scheme", y = "Accuracy (PC for continuous/count; PCCC for categorical/binary)") + theme_paper +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))
  ggsave(sprintf("figures/Fig_%s_metric1.png", sc), p, width = 7.4, height = 7.0, dpi = 300)
}
BD[, rlab := factor(ROLE_LAB[role], rev(ROLE_LAB))][, scheme := factor(scheme, SCH)][, dataset := factor(dataset, DS)]
mlab <- c(M1 = "Accuracy difference (PC or PCCC; >0 favours BMTMEMix)", M2 = "Error difference (NRMSE or Brier; <0 favours BMTMEMix)")
for (k in c("M1", "M2")) {
  p <- ggplot(BD[metric == k], aes(diff, rlab, colour = scenario)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.15, linewidth = 0.35) +
    facet_grid(dataset ~ scheme, scales = "free_x") + scale_colour_manual(values = c(S1 = "#86b6ef", S2 = "#1c5cab"), labels = c("S1 (nominal)", "S2 (ordinal)"), name = NULL) +
    labs(x = mlab[[k]], y = NULL) + theme_paper + theme(panel.grid.major.y = element_blank(), axis.text.x = element_text(size = 6))
  ggsave(sprintf("figures/Fig_diff_%s.png", k), p, width = 7.6, height = 7.6, dpi = 300)
}
if (nrow(SC)) {
SC[, scheme := factor(scheme, SCH)][, dataset := factor(dataset, DS)]
SC[, mlab2 := factor(c(M1 = "PCCC (>0 favours ordinal)", M2 = "Brier (<0 favours ordinal)", M3 = "RPS (<0 favours ordinal)")[metric],
                     c("PCCC (>0 favours ordinal)", "Brier (<0 favours ordinal)", "RPS (<0 favours ordinal)"))]
p <- ggplot(SC, aes(diff, dataset, colour = model)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.15, linewidth = 0.35) +
  facet_grid(scheme ~ mlab2, scales = "free_x") + scale_colour_manual(values = COL, name = NULL) +
  labs(x = "Ordinal minus nominal (95% cluster-bootstrap interval over lines)", y = NULL) + theme_paper + theme(axis.text.x = element_text(size = 6))
ggsave("figures/Fig_scenarios_categorical.png", p, width = 7.6, height = 6.4, dpi = 300)
}
cat("summary done\n")
