# =============================================================================
# run_all.R — reproduce every analysis of the BMTMEMix manuscript (version 3) from R.
# Portable (Windows, macOS, Linux): runs the jobs of jobs.txt sequentially, or in parallel with
#   Rscript run_all.R 2        (number of parallel workers; uses the 'parallel' package, forked on Unix,
#                               PSOCK cluster on Windows)
# The Linux/macOS helper run_all.sh does the same with xargs.
# Requirements: R >= 4.1, a C++ toolchain (Rtools on Windows, Xcode CLT on macOS), packages
#   Rcpp, RcppArmadillo, data.table, ggplot2, posterior, MASS, nnet, MCMCglmm, Matrix.
# =============================================================================
args <- commandArgs(TRUE); ncores <- if (length(args)) as.integer(args[1]) else 1L
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
dir.create("logs", showWarnings = FALSE); dir.create("outputs", showWarnings = FALSE)
rs <- file.path(R.home("bin"), "Rscript")
run <- function(cmd, log) {
  st <- system2(rs, cmd, stdout = log, stderr = log)
  cat(sprintf("%s exit=%d %s\n", paste(cmd, collapse = " "), st, format(Sys.time(), "%T")),
      file = "logs/progress.txt", append = TRUE)
  st
}
if (!file.exists("data/prepared/EYT_4.rds")) run("R/01_prepare_data.R", "logs/prepare_data.log")
run("R/01b_check_responses.R", "logs/check_responses.log")
script <- c(CV = "R/03_run_cv.R", GAUSS = "R/03_run_cv.R", SEP = "R/03_run_cv.R", REF = "R/03_run_cv.R",
            UNE = "R/03_run_cv.R", TA = "R/07_trait_assisted.R", TAX = "R/07_trait_assisted.R", TAC = "R/07b_ta_competitor.R",
            STRUCT = "R/04_structures.R", FULL = "R/05_full_fit.R", FULLLONG = "R/05b_full_fit_long.R", SENS = "R/08_sensitivity.R",
            LONG = "R/10_long_chain.R", MCMCGLMM = "R/09_mcmcglmm.R", SIM = "R/02_simulation_check.R", SBC = "R/02b_sbc.R")
jobs <- strsplit(trimws(readLines("jobs.txt")), "[[:space:]]+"); jobs <- jobs[lengths(jobs) > 0]
job_cmd <- function(j) {
  k <- j[1]; a <- j[-1]
  extra <- switch(k, GAUSS = c(a[1], a[2], "S1", "gauss"), SEP = c(a[1], "TLTE", "S1", "sep"),
                  REF = c(a[1], "TLTE", "S1", "ref"), UNE = c(a[1], a[2], "S1", "unE"), TAX = c(a[1], a[2], "xenv"), a)
  list(cmd = c(script[[k]], extra), log = file.path("logs", paste0(tolower(k), "_", paste(a, collapse = "_"), ".log")))
}
cmds <- lapply(jobs, job_cmd)
if (ncores > 1) {
  library(parallel)
  if (.Platform$OS.type == "unix") invisible(mclapply(cmds, function(x) run(x$cmd, x$log), mc.cores = ncores, mc.preschedule = FALSE))
  else { cl <- makeCluster(ncores); clusterExport(cl, c("run", "rs")); invisible(parLapplyLB(cl, cmds, function(x) run(x$cmd, x$log))); stopCluster(cl) }
} else for (x in cmds) run(x$cmd, x$log)
run("R/06_summarize.R", "logs/summarize.log")
cat("ALLDONE", format(Sys.time(), "%T"), "\n", file = "logs/progress.txt", append = TRUE)
