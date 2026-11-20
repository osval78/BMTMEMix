# =============================================================================
# 14_info_availability.R — gain of the multi-trait model as a function of the available information (v4)
# Usage: Rscript R/14_info_availability.R <dataset> <pair>
#  Bivariate BMTMEMix (target + one auxiliary trait) versus the univariate model of the target, on the ULTE line
#  folds of the main analysis (FOLD_SEED). Two factors are crossed:
#   pcand : proportion of the test (candidate) lines whose auxiliary trait is observed in every environment
#           (0 = all traits of the candidates masked, i.e. coincident masking; 1 = only the target masked);
#           nested random subsets of the candidates (25% within 50% within 100%).
#   targ  : availability of target records in training: "high" = target observed for every training line;
#           "low" = target observed for a random 25% of the training lines (the auxiliary trait stays observed
#           for all training lines, i.e. target-specific missingness in training).
#  Pairs: wheat  A = DTHD count | DTMT ordinal classes;  B = DTHD count | DTMT original (continuous);
#                C = GY continuous | DTMT original (continuous)
#         rice   A = GY continuous | GC binary;  B = GY continuous | GC original (continuous);
#                C = PHR ordinal | GC original (continuous)
#  Every prediction is saved with an indicator of whether the auxiliary trait of the candidate was observed.
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R"); source("R/v4_helpers.R")
args <- commandArgs(TRUE); d <- args[1]; pair <- args[2]
D <- readRDS(sprintf("data/prepared/%s.rds", d)); envs <- D$envOrder
eg <- eigen(D$G, symmetric = TRUE)
wheat <- grepl("^EYT", d)
spec <- if (wheat) list(A = list(tgt = "DTHD", ttype = "count", aux = "DTMT", atype = "ordinal", araw = FALSE),
                        B = list(tgt = "DTHD", ttype = "count", aux = "DTMT", atype = "continuous", araw = TRUE),
                        C = list(tgt = "GY", ttype = "continuous", aux = "DTMT", atype = "continuous", araw = TRUE)) else
                   list(A = list(tgt = "GY", ttype = "continuous", aux = "GC", atype = "binary", araw = FALSE),
                        B = list(tgt = "GY", ttype = "continuous", aux = "GC", atype = "continuous", araw = TRUE),
                        C = list(tgt = "PHR", ttype = "ordinal", aux = "GC", atype = "continuous", araw = TRUE))
sp <- spec[[pair]]
folds <- make_folds("ULTE", D$line, D$env, NFOLD, FOLD_SEED)
PC <- c(0, 0.25, 0.5, 1); TG <- c(high = 1, low = 0.25)
OUT <- list()
for (f in seq_along(folds)) {
  te <- folds[[f]]$test; train <- setdiff(seq_len(nrow(D$raw)), te)
  Yd <- make_responses(D, train)
  y_a <- if (sp$araw) D$raw[[sp$aux]] else Yd[[sp$aux]]
  set.seed(FOLD_SEED + 700 + f)
  tl <- sample(unique(D$line[te])); trl <- sample(unique(D$line[train]))
  for (tg in names(TG)) {
    keepT <- trl[seq_len(ceiling(TG[[tg]] * length(trl)))]           # training lines that keep the target
    y_t <- make_responses(D, train[D$line[train] %in% keepT])[[sp$tgt]]   # class thresholds from the retained records
    yt <- y_t; yt[te] <- NA; yt[train][!D$line[train] %in% keepT] <- NA
    # univariate model of the target (does not depend on pcand)
    set.seed(FIT_SEED + 1000 * f + 1 + (tg == "low") * 50)
    fu <- BMTMEMix(data.frame(y = yt), sp$ttype, D$line, D$env, D$G, nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
                   thin = CV_MCMC$thin, eigG = eg, envLevels = envs)
    pbu <- if (sp$ttype == "ordinal") fu$probs[[1]][fu$cell, , drop = FALSE][te, , drop = FALSE] else NULL
    for (pc in PC) {
      withAux <- tl[seq_len(round(pc * length(tl)))]
      ya <- y_a; ya[te][!D$line[te] %in% withAux] <- NA
      set.seed(FIT_SEED + 1000 * f + 10 + (tg == "low") * 50 + round(100 * pc))
      fit <- BMTMEMix(data.frame(y = yt, a = ya), c(sp$ttype, sp$atype), D$line, D$env, D$G, nIter = CV_MCMC$nIter,
                      burnIn = CV_MCMC$burnIn, thin = CV_MCMC$thin, eigG = eg, envLevels = envs)
      pbm <- if (sp$ttype == "ordinal") fit$probs[[1]][fit$cell, , drop = FALSE][te, , drop = FALSE] else NULL
      base <- data.table(dataset = d, pair = pair, target = sp$tgt, aux = paste0(sp$aux, if (sp$araw) " (original)" else paste0(" (", sp$atype, ")")),
                         type = sp$ttype, fold = f, targ = tg, pcand = pc, line = D$line[te], env = D$env[te],
                         aux_obs = D$line[te] %in% withAux, obs = y_t[te])
      m1 <- copy(base)[, `:=`(model = "MT", pred = fit$pred_sys[fit$cell, 1][te], p1 = NA_real_, p2 = NA_real_, p3 = NA_real_)]
      if (!is.null(pbm)) m1[, `:=`(pred = NA_real_, p1 = pbm[, 1], p2 = pbm[, 2], p3 = pbm[, 3])]
      u1 <- copy(base)[, `:=`(model = "UNI", pred = fu$pred_sys[fu$cell, 1][te], p1 = NA_real_, p2 = NA_real_, p3 = NA_real_)]
      if (!is.null(pbu)) u1[, `:=`(pred = NA_real_, p1 = pbu[, 1], p2 = pbu[, 2], p3 = pbu[, 3])]
      gt <- copy(base)[, `:=`(model = "MT genetic value", pred = fit$gv[fit$cell, 1][te], p1 = NA_real_, p2 = NA_real_, p3 = NA_real_)]
      OUT[[length(OUT) + 1]] <- rbind(m1, u1, if (sp$ttype != "ordinal") gt)
    }
  }
  cat("fold", f, "done", format(Sys.time()), "\n")
}
dir.create(file.path(V4DIR, "info"), showWarnings = FALSE, recursive = TRUE)
saveRDS(rbindlist(OUT), file.path(V4DIR, "info", sprintf("%s_%s.rds", d, pair)))
cat("done", d, pair, "\n")
