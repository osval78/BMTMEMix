# Usage: Rscript R/01b_check_responses.R
# Checks that the responses rebuilt by make_responses() with all records as training set (the function used in
# cross-validation with training-set thresholds) are identical to the responses stored in data/prepared/*.rds.
source("R/00_settings.R")
make_responses <- local({ src <- readLines("R/BMTMEMix.R"); i <- grep("^make_responses <- function", src)
  j <- i + which(src[(i + 1):length(src)] == "}")[1]; eval(parse(text = src[i:j])); make_responses })
out <- do.call(rbind, lapply(DATASETS, function(d) {
  D <- readRDS(sprintf("data/prepared/%s.rds", d)); Y <- make_responses(D)
  data.frame(dataset = d, identical = isTRUE(all.equal(as.data.frame(Y), as.data.frame(D$Yd), check.attributes = FALSE)),
             n = nrow(Y)) }))
print(out); dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "outputs/tables/check_responses.csv", row.names = FALSE)
