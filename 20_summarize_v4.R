# =============================================================================
# 20_summarize_v4.R — tables (outputs/tables/v4_*.csv) and figures of the version-4 analyses
# Requires the outputs of R/11–R/19 (outputs/v4/...) and of R/06_summarize.R (cv_predictions_all.csv.gz).
# Every accuracy is computed within environment on pooled out-of-fold predictions and averaged over environments,
# and every difference is accompanied by a paired cluster-bootstrap interval over lines (BMTME_BOOT replicates).
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(parallel) })
source("R/00_settings.R"); source("R/v4_helpers.R")
want <- function(x) { s <- Sys.getenv("BMTME_SECTIONS", "all"); s == "all" || x %in% strsplit(s, ",")[[1]] }
TB <- "outputs/tables"; NC <- as.integer(Sys.getenv("BMTME_CORES", "2"))
theme_paper <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
        strip.text = element_text(face = "bold"), legend.position = "top",
        plot.background = element_rect(fill = "white", colour = NA))
rdall <- function(dir) { f <- list.files(file.path(V4DIR, dir), pattern = "\\.rds$", full.names = TRUE); if (!length(f)) return(NULL)
  rbindlist(lapply(f, readRDS), fill = TRUE) }
pooled <- function(dt, predcol = "pred") {     # metrics of one model on pooled records (within env, averaged)
  rl <- role_of(dt$type[1])
  ee <- if ("env" %in% names(dt)) dt$env else rep(1L, nrow(dt))
  v <- sapply(split(seq_len(nrow(dt)), ee), function(ix) { x <- dt[ix]; mvals(x$obs, x[[predcol]], cbind(x$p1, x$p2, x$p3), rl) })
  r <- rowMeans(matrix(v, nrow = 3), na.rm = TRUE); list(M1 = r[1], M2 = r[2], M3 = r[3])
}
bootjobs <- function(jobs) rbindlist(mclapply(jobs, function(j) { b <- boot_pair(j$w); for (k in names(j$key)) b[[k]] <- j$key[[k]]; b }, mc.cores = NC), fill = TRUE)

## ============================ 1. Monte Carlo precision and matched priors ============================
MC <- rdall("mc_prior")
if (!is.null(MC) && want("mc")) {
  MC[, role := role_of(type)]
  MC[, tgt := ifelse(is.na(target), trait, target)]
  MC <- MC[is.na(target) | trait == target]
  # pooled accuracy by stream (systematic predictor for continuous/count; probabilities for discrete traits)
  MC[, predsel := ifelse(role %in% c("continuous", "count"), pred_sys, pred)]
  acc <- MC[, pooled(.SD, "predsel"), by = .(dataset, scenario, scheme, trait, role, model, stream)]
  fwrite(acc, file.path(TB, "v4_mc_accuracy_by_stream.csv"))
  w <- dcast(acc, dataset + scenario + scheme + trait + role + stream ~ model, value.var = "M1")
  w[, `:=`(d_MT = MT - UNI, d_MTpm = `MT-pm` - UNI, d_pm = `MT-pm` - MT)]
  mcs <- w[, .(UNI = mean(UNI), MT = mean(MT), MTpm = mean(`MT-pm`), sd_UNI = sd(UNI), sd_MT = sd(MT),
               d_MT = mean(d_MT), sd_d_MT = sd(d_MT), d_MTpm = mean(d_MTpm), sd_d_MTpm = sd(d_MTpm), d_pm = mean(d_pm),
               range_d_MT = max(d_MT) - min(d_MT)), by = .(dataset, scenario, scheme, trait, role)]
  rel <- MC[, .(rel_mcse_median = median(mcse_eta / sd_eta, na.rm = TRUE), rel_mcse_p95 = quantile(mcse_eta / sd_eta, 0.95, na.rm = TRUE),
                ess_median = median(ess_eta, na.rm = TRUE)), by = .(dataset, scenario, scheme, trait, model)]
  relw <- dcast(rel, dataset + scenario + scheme + trait ~ model, value.var = "rel_mcse_median")
  setnames(relw, c("MT", "MT-pm", "UNI"), c("relMCSE_MT", "relMCSE_MTpm", "relMCSE_UNI"), skip_absent = TRUE)
  mcs <- merge(mcs, relw, by = c("dataset", "scenario", "scheme", "trait"), all.x = TRUE)
  # bootstrap of matched-prior MT minus UNI (stream 1) and of MT minus UNI (stream 1)
  jobs <- list()
  for (k in unique(MC[, paste(dataset, scenario, scheme, trait, sep = "|")])) {
    z <- strsplit(k, "\\|")[[1]]; dt <- MC[dataset == z[1] & scenario == z[2] & scheme == z[3] & trait == z[4] & stream == 1]
    dt[, pred := predsel]
    for (m in c("MT", "MT-pm")) jobs[[length(jobs) + 1]] <- list(w = widen2(dt, "UNI", m), key = list(dataset = z[1], scenario = z[2], scheme = z[3], trait = z[4], model = m))
  }
  BB <- bootjobs(jobs)[metric == "M1"]
  bw <- dcast(BB, dataset + scenario + scheme + trait ~ model, value.var = c("lo", "hi"))
  mcs <- merge(mcs, bw, by = c("dataset", "scenario", "scheme", "trait"), all.x = TRUE)
  fwrite(mcs, file.path(TB, "v4_mc_prior_summary.csv"))
  fwrite(rel, file.path(TB, "v4_mc_relative_mcse.csv"))
}

## ============================ 2. repeated partitions ============================
PT <- list.files(file.path(V4DIR, "partitions"), pattern = "\\.rds$", full.names = TRUE)
P0 <- if (file.exists(file.path(TB, "cv_predictions_all.csv.gz"))) fread(file.path(TB, "cv_predictions_all.csv.gz")) else NULL
if (length(PT) && !is.null(P0) && want("part")) {
  PP <- rbindlist(lapply(PT, function(f) as.data.table(readRDS(f)$preds)), fill = TRUE)
  P0s <- P0[scenario == "S2" & scheme %in% c("TLTE", "ULTE")][, partition := 2026L]
  PP <- rbind(PP, P0s[dataset %in% unique(PP$dataset)], fill = TRUE)
  PP[, role := role_of(type)]
  acc <- PP[, pooled(.SD), by = .(dataset, scheme, trait, role, model, partition)]
  w <- dcast(acc, dataset + scheme + trait + role + partition ~ model, value.var = "M1")[, d := BMTMEMix - Univariate]
  fwrite(w, file.path(TB, "v4_partitions_by_seed.csv"))
  fwrite(w[, .(UNI = mean(Univariate), MT = mean(BMTMEMix), d_mean = mean(d), d_sd = sd(d), d_min = min(d), d_max = max(d), n = .N),
           by = .(dataset, scheme, trait, role)], file.path(TB, "v4_partitions_summary.csv"))
}

## ============================ 3. trait-assisted with paired intervals ============================
TA <- rdall("ta")
if (!is.null(TA) && want("ta")) {
  TA[, role := role_of(type)]
  # univariate model and competitors of S2 for non-categorical targets are those of S1 (identical data and seeds)
  reuse <- TA[scenario == "S1" & scheme == "ULTE-TA" & model %in% c("UNI", "UNI+sec (same plot)", "UNI+sec (across env, ridge)") & role != "categorical"]
  TA <- rbind(TA, copy(reuse)[, scenario := "S2"])
  acc <- TA[, pooled(.SD), by = .(dataset, scenario, scheme, trait, role, model)]
  fwrite(acc, file.path(TB, "v4_ta_accuracy.csv"))
  cmp <- list(c("UNI", "MT-TA systematic"), c("UNI", "MT-TA phenotypic"), c("UNI+sec (same plot)", "MT-TA systematic"),
              c("UNI+sec (across env, ridge)", "MT-TA systematic"), c("UNI", "UNI+sec (same plot)"), c("UNI", "UNI+sec (across env, ridge)"),
              c("MT-TA systematic", "MT-TA phenotypic"))
  jobs <- list()
  for (k in unique(TA[scheme == "ULTE-TA", paste(dataset, scenario, trait, sep = "|")])) {
    z <- strsplit(k, "\\|")[[1]]; dt <- TA[scheme == "ULTE-TA" & dataset == z[1] & scenario == z[2] & trait == z[3]]
    for (cc in cmp) if (all(cc %in% dt$model)) jobs[[length(jobs) + 1]] <- list(w = widen2(dt, cc[1], cc[2]), key = list(dataset = z[1], scenario = z[2], scheme = "ULTE-TA", trait = z[3], mA = cc[1], mB = cc[2]))
  }
  # xenv: systematic predictor versus UNI on the same records
  for (k in unique(TA[scheme == "ULTE-TA-xenv", paste(dataset, scenario, trait, sep = "|")])) {
    z <- strsplit(k, "\\|")[[1]]
    dt <- rbind(TA[scheme == "ULTE-TA-xenv" & dataset == z[1] & scenario == z[2] & trait == z[3] & model == "MT-TA systematic"][, model := "xenv systematic"],
                TA[scheme == "ULTE-TA" & dataset == z[1] & scenario == z[2] & trait == z[3] & model %in% c("UNI", "MT-TA systematic")], fill = TRUE)
    jobs[[length(jobs) + 1]] <- list(w = widen2(dt, "UNI", "xenv systematic"), key = list(dataset = z[1], scenario = z[2], scheme = "ULTE-TA-xenv", trait = z[3], mA = "UNI", mB = "xenv systematic"))
    jobs[[length(jobs) + 1]] <- list(w = widen2(dt, "xenv systematic", "MT-TA systematic"), key = list(dataset = z[1], scenario = z[2], scheme = "ULTE-TA-xenv", trait = z[3], mA = "xenv systematic", mB = "MT-TA systematic"))
  }
  TB4 <- bootjobs(jobs)
  TB4 <- merge(TB4, unique(TA[, .(dataset, trait, role)]), by = c("dataset", "trait"))
  TB4 <- TB4[!(role %in% c("continuous", "count", "binary") & metric == "M3")]
  fwrite(TB4, file.path(TB, "v4_ta_bootstrap.csv"))
  fwrite(TB4[, .(mean_diff = mean(diff), min_diff = min(diff), max_diff = max(diff), n = .N,
                 n_pos = sum((metric == "M1" & lo > 0) | (metric != "M1" & hi < 0)), n_neg = sum((metric == "M1" & hi < 0) | (metric != "M1" & lo > 0))),
             by = .(scheme, mA, mB, role, metric)], file.path(TB, "v4_ta_bootstrap_summary.csv"))
  # figure: MT-TA systematic minus UNI and minus the two competitors (accuracy)
  fd <- TB4[metric == "M1" & mB == "MT-TA systematic" & scheme == "ULTE-TA"]
  fd[, cmpl := factor(c("UNI" = "vs UNI", "UNI+sec (same plot)" = "vs same-plot regression", "UNI+sec (across env, ridge)" = "vs across-env. ridge")[mA],
                      c("vs UNI", "vs same-plot regression", "vs across-env. ridge"))]
  fd[, rlab := factor(c(continuous = "Continuous", count = "Count", categorical = "Categorical", binary = "Binary")[role], rev(c("Continuous", "Count", "Categorical", "Binary")))]
  fd[, dataset := factor(dataset, DATASETS)]
  p <- ggplot(fd, aes(diff, rlab, colour = scenario)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.15, linewidth = 0.35) +
    facet_grid(dataset ~ cmpl, scales = "free_x") + scale_colour_manual(values = c(S1 = "#86b6ef", S2 = "#1c5cab"), labels = c("S1 (nominal)", "S2 (ordinal)"), name = NULL) +
    labs(x = "MT-TA (systematic) minus comparator, accuracy (PC or PCCC), 95% paired interval", y = NULL) +
    theme_paper + theme(panel.grid.major.y = element_blank(), axis.text.x = element_text(size = 6))
  ggsave("figures/Fig_TA_paired.png", p, width = 7.6, height = 7.4, dpi = 300)
}

## ============================ 4. availability of information ============================
IN <- rdall("info")
if (!is.null(IN) && want("info")) {
  IN[, role := role_of(type)]
  sub <- rbind(copy(IN)[, grp := "all candidates"], IN[pcand > 0 & pcand < 1 & aux_obs == TRUE][, grp := "with auxiliary"],
               IN[pcand > 0 & pcand < 1 & aux_obs == FALSE][, grp := "without auxiliary"])
  acc <- sub[, pooled(.SD), by = .(dataset, pair, target, aux, role, targ, pcand, grp, model)]
  fwrite(acc, file.path(TB, "v4_info_accuracy.csv"))
  jobs <- list()
  for (k in unique(sub[, paste(dataset, pair, targ, pcand, grp, sep = "|")])) {
    z <- strsplit(k, "\\|")[[1]]; dt <- sub[dataset == z[1] & pair == z[2] & targ == z[3] & pcand == as.numeric(z[4]) & grp == z[5]]
    jobs[[length(jobs) + 1]] <- list(w = widen2(dt, "UNI", "MT"), key = list(dataset = z[1], pair = z[2], targ = z[3], pcand = as.numeric(z[4]), grp = z[5]))
  }
  IB <- bootjobs(jobs)
  IB <- merge(IB, unique(IN[, .(dataset, pair, target, aux, role)]), by = c("dataset", "pair"))[metric == "M1" | (role == "categorical" & metric == "M3") | metric == "M2"]
  fwrite(IB, file.path(TB, "v4_info_bootstrap.csv"))
  # per environment (MT minus UNI accuracy), all candidates
  pe <- IN[, pooled(.SD), by = .(dataset, pair, target, aux, role, targ, pcand, env, model)]
  pew <- dcast(pe, dataset + pair + target + aux + role + targ + pcand + env ~ model, value.var = "M1")[, d := MT - UNI]
  fwrite(pew, file.path(TB, "v4_info_by_environment.csv"))
  f <- IB[metric == "M1" & grp == "all candidates"]
  f[, lab := paste0(dataset, ": ", target, " | ", aux)]
  f[, targ := factor(targ, c("high", "low"), c("target in all training lines", "target in 25% of training lines"))]
  p <- ggplot(f, aes(100 * pcand, diff, colour = targ)) + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.5) + geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.2, linewidth = 0.35, position = position_dodge(width = 5)) +
    facet_wrap(~lab, ncol = 3, scales = "free_y") + scale_colour_manual(values = c("#1c5cab", "#eb6834"), name = NULL) +
    scale_x_continuous(breaks = c(0, 25, 50, 100)) +
    labs(x = "Candidates with the auxiliary trait observed (%)", y = "MT minus UNI accuracy (PC or PCCC), 95% paired interval") + theme_paper +
    theme(strip.text = element_text(size = 7))
  ggsave("figures/Fig_info_availability.png", p, width = 7.6, height = 6.6, dpi = 300)
}

## ============================ 5. Gaussian approximation, comparable uncertainty ============================
GA <- rdall("gauss")
if (!is.null(GA) && !is.null(P0) && want("gauss")) {
  GA[, role := role_of(type)]
  B0 <- P0[model == "BMTMEMix" & scheme %in% c("TLTE", "ULTE") & role %in% c("categorical", "binary"),
           .(dataset, scheme, fold, trait, type, line, env, obs, pred, p1, p2, p3, role, model = paste0("BMTMEMix ", scenario))]
  AL <- rbind(GA, B0, fill = TRUE)
  acc <- AL[, pooled(.SD), by = .(dataset, scheme, trait, role, model)]
  fwrite(acc, file.path(TB, "v4_gauss_accuracy.csv"))
  jobs <- list()
  for (k in unique(AL[, paste(dataset, scheme, trait, sep = "|")])) {
    z <- strsplit(k, "\\|")[[1]]; dt <- AL[dataset == z[1] & scheme == z[2] & trait == z[3]]
    for (bm in intersect(c("BMTMEMix S1", "BMTMEMix S2"), unique(dt$model))) for (gm in c("Gaussian (plug-in)", "Gaussian (draw-wise)", "Gaussian (draw-wise, training-calibrated)"))
      jobs[[length(jobs) + 1]] <- list(w = widen2(dt, gm, bm), key = list(dataset = z[1], scheme = z[2], trait = z[3], mA = gm, mB = bm))
  }
  GB <- bootjobs(jobs); GB <- merge(GB, unique(AL[, .(dataset, trait, role)]), by = c("dataset", "trait"))
  GB <- GB[!(role == "binary" & metric == "M3")]
  fwrite(GB, file.path(TB, "v4_gauss_bootstrap.csv"))
  fwrite(GB[, .(gauss = mean(A), mix = mean(B), mean_diff = mean(diff), n = .N,
                n_mix_better = sum((metric == "M1" & lo > 0) | (metric != "M1" & hi < 0)), n_gauss_better = sum((metric == "M1" & hi < 0) | (metric != "M1" & lo > 0))),
             by = .(role, mA, mB, metric)], file.path(TB, "v4_gauss_summary.csv"))
}

## ============================ 6. untested environments ============================
UN <- rdall("untested")
if (!is.null(UN) && want("unt")) {
  UN[, role := role_of(type)]
  ue <- UN[role %in% c("continuous", "count"), .(PC = suppressWarnings(cor(obs, pred)), bias = mean(pred - obs), rel_bias = mean(pred - obs) / mean(obs),
              MAE = mean(abs(pred - obs)), slope = if (sd(pred) > 0) coef(lm(obs ~ pred))[2] else NA_real_,
              cover_q = mean(obs >= q5 & obs <= q95), cover_sd = mean(abs(obs - pred) <= 1.645 * pred_sd), n = .N),
           by = .(dataset, scheme, env, trait, role, model)]
  fwrite(ue, file.path(TB, "v4_untested_by_environment.csv"))
  fwrite(ue[, .(PC = mean(PC), abs_rel_bias = mean(abs(rel_bias)), MAE_rel = NA_real_, slope = median(slope, na.rm = TRUE),
                cover_q = mean(cover_q), cover_sd = mean(cover_sd), n_env = .N), by = .(scheme, role, model)], file.path(TB, "v4_untested_summary.csv"))
  jobs <- list()
  for (k in unique(UN[, paste(dataset, scheme, trait, sep = "|")])) {
    z <- strsplit(k, "\\|")[[1]]; dt <- UN[dataset == z[1] & scheme == z[2] & trait == z[3]]
    for (m in c("MT", "MT-GxEdiag")) jobs[[length(jobs) + 1]] <- list(w = widen2(dt, "UNI", m), key = list(dataset = z[1], scheme = z[2], trait = z[3], model = m))
    jobs[[length(jobs) + 1]] <- list(w = widen2(dt, "MT", "MT-GxEdiag"), key = list(dataset = z[1], scheme = z[2], trait = z[3], model = "GxEdiag minus MT"))
  }
  UB <- bootjobs(jobs); UB <- merge(UB, unique(UN[, .(dataset, trait, role)]), by = c("dataset", "trait"))[metric %in% c("M1", "M2")]
  fwrite(UB, file.path(TB, "v4_untested_bootstrap.csv"))
}

## ============================ 7. simulation: ideal vs estimated ============================
OR <- { f <- list.files(file.path(V4DIR, "sim_oracle"), full.names = TRUE); if (length(f)) rbindlist(lapply(f, fread)) else NULL }
if (!is.null(OR) && want("oracle")) {
  fwrite(OR[, .(acc = mean(acc), se = sd(acc) / sqrt(.N), n = .N), by = .(ttype, h2, rho, mask, model)], file.path(TB, "v4_oracle_accuracy.csv"))
  w <- dcast(OR, cfg + ttype + h2 + rho + rep + mask ~ model, value.var = "acc")
  cont <- list(c("MT-true", "UNI-true", "ideal: MT - UNI"), c("MT-est", "UNI-est", "estimated: MT - UNI"),
               c("MT-est", "MT-true", "cost of estimation, MT"), c("UNI-est", "UNI-true", "cost of estimation, UNI"))
  OD <- rbindlist(lapply(cont, function(cc) w[, .(contrast = cc[3], diff = mean(get(cc[1]) - get(cc[2])), se = sd(get(cc[1]) - get(cc[2])) / sqrt(.N), n = .N),
                                              by = .(ttype, h2, rho, mask)]))
  fwrite(OD, file.path(TB, "v4_oracle_contrasts.csv"))
  OD[, lab := factor(paste0(ttype, " target, h2 = ", h2))]
  p <- ggplot(OD[contrast %in% c("ideal: MT - UNI", "estimated: MT - UNI")], aes(rho, diff, colour = contrast, linetype = mask)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) + geom_line(linewidth = 0.5) +
    geom_pointrange(aes(ymin = diff - 1.96 * se, ymax = diff + 1.96 * se), size = 0.2, linewidth = 0.35) +
    facet_wrap(~lab, ncol = 4) + scale_colour_manual(values = c("ideal: MT - UNI" = "#1c5cab", "estimated: MT - UNI" = "#eb6834"), name = NULL) + scale_x_continuous(breaks = c(0, 0.5, 0.8)) +
    labs(x = "Genetic correlation between target and secondary trait", y = "MT - UNI accuracy of genetic values", linetype = NULL) +
    theme_paper + theme(legend.box = "vertical")
  ggsave("figures/Fig_oracle.png", p, width = 7.6, height = 4.2, dpi = 300)
}

## ============================ 8. simulation: changing correlation in a new environment ============================
EC <- { f <- list.files(file.path(V4DIR, "sim_envcor"), full.names = TRUE); if (length(f)) rbindlist(lapply(f, fread)) else NULL }
if (!is.null(EC) && want("envcor")) fwrite(EC[, .(PC_pheno = mean(PC_pheno), se_PC = sd(PC_pheno) / sqrt(.N), PC_g = mean(PC_g), bias = mean(bias), abs_bias = mean(abs(bias)),
                                MAE = mean(MAE), n = .N), by = .(scenario, mask, model)], file.path(TB, "v4_envcor_summary.csv"))

## ============================ 9. order of traits and environments ============================
OO <- { f <- list.files(file.path(V4DIR, "order"), full.names = TRUE); if (length(f)) rbindlist(lapply(f, fread), fill = TRUE) else NULL }
if (!is.null(OO) && want("order")) {
  om <- OO[, .(M1 = mean(M1)), by = .(dataset, analysis, structure, order, justified, original)]
  fwrite(om, file.path(TB, "v4_order_by_order.csv"))
  fwrite(om[, .(n_orders = .N, min = min(M1), max = max(M1), range = max(M1) - min(M1), sd = sd(M1), justified = M1[justified][1],
                original = mean(M1[original])), by = .(dataset, analysis, structure)], file.path(TB, "v4_order_summary.csv"))
  ot <- OO[, .(range = max(M1) - min(M1)), by = .(dataset, analysis, structure, trait)]
  fwrite(ot, file.path(TB, "v4_order_range_by_trait.csv"))
}

## ============================ 10. selection decisions: top 20% for grain yield ============================
topagree <- function(o, p, q = 0.2) { k <- max(1, round(q * length(o))); mean(order(-p)[1:k] %in% order(-o)[1:k]) }
if (!is.null(P0) && want("sel")) {
  gy <- P0[role == "continuous" & scheme %in% c(SCHEMES, "CHRONO")]
  sel <- gy[, .(top20 = topagree(obs, pred), n = .N), by = .(dataset, scenario, scheme, model, env)][, .(top20 = mean(top20), n_env = .N), by = .(dataset, scenario, scheme, model)]
  fwrite(sel, file.path(TB, "v4_selection_top20_cv.csv"))
}
if (!is.null(TA) && want("sel")) {
  gy <- TA[role == "continuous" & scheme == "ULTE-TA"]
  sel <- gy[, .(top20 = topagree(obs, pred), n = .N), by = .(dataset, scenario, model, env)][, .(top20 = mean(top20), n_env = .N), by = .(dataset, scenario, model)]
  fwrite(sel, file.path(TB, "v4_selection_top20_ta.csv"))
}
## ============================ 11. additional SBC ============================
SB <- { f <- list.files(file.path(V4DIR, "sbc"), full.names = TRUE); if (length(f)) rbindlist(lapply(f, fread)) else NULL }
if (!is.null(SB) && want("sbc")) {
  chi <- SB[, { tab <- table(cut(rank, breaks = seq(-0.5, 99.5, length.out = 11))); ct <- suppressWarnings(chisq.test(tab)); .(chisq = unname(ct$statistic), p_value = ct$p.value, n = .N) },
            by = .(config, px, parameter)]
  fwrite(chi, file.path(TB, "v4_sbc_extra.csv"))
  hl <- SB[, .(h = .N / 10), by = .(config, px, parameter)]
  SB[, lab := paste0(config, ifelse(px, " + scale move", ""))]; hl[, lab := paste0(config, ifelse(px, " + scale move", ""))]
  p <- ggplot(SB, aes(rank)) + geom_histogram(breaks = seq(-0.5, 99.5, length.out = 11), fill = "#2a78d6", colour = "white") +
    geom_hline(data = hl, aes(yintercept = h), colour = "#eb6834", linewidth = 0.4) + facet_wrap(lab ~ parameter, ncol = 9, scales = "free_y") +
    labs(x = "Rank of the true value among 99 posterior draws", y = "Replicates") + theme_paper + theme(strip.text = element_text(size = 5), axis.text = element_text(size = 5))
  ggsave("figures/FigS_sbc_extra.png", p, width = 10, height = 4.2, dpi = 200)
}
cat("v4 summary done\n")
