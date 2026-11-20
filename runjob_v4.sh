#!/bin/bash
# runjob_v4.sh <JOB> <args...> — one version-4 analysis job (see jobs_v4.txt)
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p logs/v4
k=$1; a=$2; b=$3; c=$4
# skip jobs already completed (allows the queue to be restarted)
if grep -q "^$k $a $b $c $5 exit=0" logs/v4/progress.txt 2>/dev/null; then exit 0; fi
case $k in
  MCP)    Rscript R/11_mc_prior.R $a $b $c          > logs/v4/mcp_${a}_${b}_${c}.log 2>&1 ;;
  PART)   Rscript R/12_partitions.R $a $b $c        > logs/v4/part_${a}_${b}_${c}.log 2>&1 ;;
  TA4)    Rscript R/13_ta_v4.R $a $b $c             > logs/v4/ta4_${a}_${b}_${c}.log 2>&1 ;;
  INFO)   Rscript R/14_info_availability.R $a $b    > logs/v4/info_${a}_${b}.log 2>&1 ;;
  GAUSS4) Rscript R/15_gauss_v4.R $a $b             > logs/v4/gauss4_${a}_${b}.log 2>&1 ;;
  UNT4)   Rscript R/16_untested_v4.R $a $b          > logs/v4/unt4_${a}_${b}.log 2>&1 ;;
  ORACLE) Rscript R/17_sim_oracle.R $a $b $c        > logs/v4/oracle_${a}_${b}.log 2>&1 ;;
  ENVCOR) Rscript R/18_sim_envcor.R $a $b           > logs/v4/envcor_${a}.log 2>&1 ;;
  ORDER)  Rscript R/19_order.R $a $b                > logs/v4/order_${a}_${b}.log 2>&1 ;;
  SBC4)   Rscript R/02c_sbc_extra.R $a $b $c $5      > logs/v4/sbc4_${a}_${b}_$5.log 2>&1 ;;
esac
echo "$k $a $b $c $5 exit=$? $(date +%T)" >> logs/v4/progress.txt
