#!/bin/bash
# Reproduce all analyses of the BMTMEMix manuscript (version 3) on Linux/macOS with two jobs in parallel.
# Portable alternative for any operating system: Rscript run_all.R 2
set -e
cd "$(dirname "$0")"
mkdir -p logs outputs
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
[ -f data/prepared/EYT_4.rds ] || Rscript R/01_prepare_data.R > logs/prepare_data.log 2>&1
Rscript R/01b_check_responses.R > logs/check_responses.log 2>&1
: > logs/progress.txt
cat jobs.txt | xargs -P 2 -L 1 ./runjob.sh
echo "ALLDONE $(date +%T)" >> logs/progress.txt
Rscript R/06_summarize.R > logs/summarize.log 2>&1
