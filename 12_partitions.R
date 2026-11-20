# =============================================================================
# 12_partitions.R — repeated random partitions (v4)
# Usage: Rscript R/12_partitions.R <dataset> <TLTE|ULTE> <partition_seed>
#  Scenario 2; BMTMEMix and the four univariate models with a new random partition (seed 2027 or 2028 instead of
#  FOLD_SEED = 2026); fit seeds as in the main analysis. Used to show how much the MT - UNI differences change from
#  one partition to another (the bootstrap intervals are conditional on one partition).
# =============================================================================
source("R/00_settings.R"); source("R/BMTMEMix.R"); source("R/v4_helpers.R")
args <- commandArgs(TRUE); d <- args[1]; sch <- args[2]; ps <- as.integer(args[3])
D <- readRDS(sprintf("data/prepared/%s.rds", d)); types <- D$types_S2; eg <- eigen(D$G, symmetric = TRUE)
out <- run_cv(D, types, sch, nfold = NFOLD, seed = ps, covT = "UN", covE = covE_for(sch), nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn,
              thin = CV_MCMC$thin, eigG = eg, fit_seed = FIT_SEED)
out$preds$dataset <- d; out$preds$partition <- ps; out$metrics$dataset <- d; out$metrics$partition <- ps
dir.create(file.path(V4DIR, "partitions"), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, file.path(V4DIR, "partitions", sprintf("%s_%s_%d.rds", d, sch, ps)))
cat("done", d, sch, ps, "\n")
