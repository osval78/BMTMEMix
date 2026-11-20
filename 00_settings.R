# Global settings shared by all analysis scripts (version 3)
DATASETS  <- c("EYT_4", "EYT_5", "EYT_6", "EYT_7", "Indica", "Japonica")
SCHEMES   <- c("TLTE", "ULTE", "TLUE", "ULUE")
NFOLD     <- 5          # folds for TLTE/ULTE and line groups for ULUE (crossed with every environment)
FOLD_SEED <- 2026       # partitions (identical for every model, scenario and variant)
FIT_SEED  <- 2026       # base seed of the model fits (see run_cv)
CV_MCMC   <- list(nIter = 2500, burnIn = 700, thin = 4)    # cross-validation and trait-assisted (450 draws)
STR_MCMC  <- list(nIter = 2500, burnIn = 700, thin = 4)    # covariance structures (same precision as CV)
FULL_MCMC <- list(nIter = 6000, burnIn = 2000, thin = 4, chains = 3, initDisp = 1)  # full fits, 3 dispersed chains
UNTESTED_COVE <- "ID"   # exchangeable G x E for untested environments (TLUE, ULUE, CHRONO)
scenario_types <- function(D, scen) if (scen == "S1") D$types else D$types_S2
covE_for <- function(scheme) if (scheme %in% c("TLUE", "ULUE", "CHRONO")) UNTESTED_COVE else "UN"
