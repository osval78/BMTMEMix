#!/bin/bash
# runjob.sh <JOB> <args...> — one analysis job (see jobs.txt)
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
k=$1; a=$2; b=$3; c=$4
case $k in
  CV)       Rscript R/03_run_cv.R $a $b $c            > logs/cv_${a}_${b}_${c}.log 2>&1 ;;
  GAUSS)    Rscript R/03_run_cv.R $a $b S1 gauss      > logs/gauss_${a}_${b}.log 2>&1 ;;
  SEP)      Rscript R/03_run_cv.R $a TLTE S1 sep      > logs/sep_${a}.log 2>&1 ;;
  REF)      Rscript R/03_run_cv.R $a TLTE S1 ref      > logs/ref_${a}.log 2>&1 ;;
  UNE)      Rscript R/03_run_cv.R $a $b S1 unE        > logs/une_${a}_${b}.log 2>&1 ;;
  TA)       Rscript R/07_trait_assisted.R $a $b       > logs/ta_${a}_${b}.log 2>&1 ;;
  TAX)      Rscript R/07_trait_assisted.R $a $b xenv  > logs/tax_${a}_${b}.log 2>&1 ;;
  TAC)      Rscript R/07b_ta_competitor.R $a $b       > logs/tac_${a}_${b}.log 2>&1 ;;
  STRUCT)   Rscript R/04_structures.R $a              > logs/struct_${a}.log 2>&1 ;;
  FULL)     Rscript R/05_full_fit.R $a $b             > logs/full_${a}_${b}.log 2>&1 ;;
  FULLLONG) Rscript R/05b_full_fit_long.R $a $b       > logs/fulllong_${a}_${b}.log 2>&1 ;;
  SENS)     Rscript R/08_sensitivity.R $a             > logs/sens_${a}.log 2>&1 ;;
  LONG)     Rscript R/10_long_chain.R $a              > logs/long_${a}.log 2>&1 ;;
  MCMCGLMM) Rscript R/09_mcmcglmm.R $a $b             > logs/mcmcglmm_${a}.log 2>&1 ;;
  SIM)      Rscript R/02_simulation_check.R $a $b $c  > logs/sim_${a}_${b}.log 2>&1 ;;
  SBC)      Rscript R/02b_sbc.R $a $b $c              > logs/sbc${c}_${a}.log 2>&1 ;;
esac
echo "$k $a $b $c exit=$? $(date +%T)" >> logs/progress.txt
