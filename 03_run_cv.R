# Usage: Rscript R/03_run_cv.R <dataset> <scheme> <S1|S2> [variant]
#   S1: continuous, count, nominal, binary -> BMTMEMix + four univariate models
#   S2: continuous, count, ordinal, binary -> BMTMEMix + univariate ordinal model (the other three
#       univariate models are identical to S1 and are not refitted)
#   variants: gauss (all traits treated as Gaussian, multi-trait; competitor)
#             sep   (purely separable genetic covariance, multi-trait, TLTE)
#             unE   (untested environments with unstructured Omega, BMTMEMix + univariate)
#             ref   (nominal model with class 1 as reference category, multi-trait, TLTE)
source("R/00_settings.R"); source("R/BMTMEMix.R")
args <- commandArgs(TRUE); d <- args[1]; sch <- args[2]; scen <- args[3]
variant <- if (length(args) > 3) args[4] else ""
D <- readRDS(sprintf("data/prepared/%s.rds", d)); types <- scenario_types(D, scen)
eg <- eigen(D$G, symmetric = TRUE)
covE <- if (variant == "unE") "UN" else covE_for(sch)
gxe <- if (variant == "sep") "separable" else "main"
uni <- !(variant %in% c("gauss", "sep", "ref"))
ut  <- if (scen == "S2") which(types == "ordinal") else NULL
relabel <- if (variant == "ref") list(trait = which(types == "nominal"), map = c(3L, 2L, 1L)) else NULL
out <- run_cv(D, types, sch, nfold = NFOLD, seed = FOLD_SEED, gxe = gxe, covT = "UN", covE = covE,
              nIter = CV_MCMC$nIter, burnIn = CV_MCMC$burnIn, thin = CV_MCMC$thin, eigG = eg,
              univariate = uni, uni_traits = ut, fit_seed = FIT_SEED, relabel = relabel, gaussian = (variant == "gauss"))
out$metrics$dataset <- d; out$metrics$scenario <- scen; out$metrics$variant <- ifelse(variant == "", "main", variant)
out$metrics$covE <- covE; out$preds$dataset <- d; out$preds$scenario <- scen
dir <- file.path("outputs", "cv", if (variant == "") scen else paste0(scen, "_", variant))
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(out, file.path(dir, sprintf("%s_%s.rds", d, sch)))
write.csv(out$metrics, file.path(dir, sprintf("%s_%s_metrics.csv", d, sch)), row.names = FALSE)
cat("done", d, sch, scen, variant, "\n")
