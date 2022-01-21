#!/bin/bash

trap "exit" INT TERM ERR
trap "kill 0" EXIT

# mkdir ~/output
# rsync /proj/snic2021-22-887/bayesian_correlations/bayesian_correlations.sif

current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
run_dir=$current_dir/run_dir
build_dir=$current_dir/build_dir

# args
outer_seed=2021
sample_out=200
initial_size=10
step_size=10
n_sim=100
stan_chains=4
stan_cores=1
stan_warmup=1000
stan_iter=5000
mapply_cores=23

# -------------------------------
# changes for each model go here!
# -------------------------------

model='weakly_inf'
pop_rho='0.3'

# -------------
# run container
# -------------

# test if this will work:
# singularity run -B ~/output:/output bayesian_correlations.sif \

start=$(date +%s)
singularity run -B $run_dir/output:/output $build_dir/bayesian_correlations.sif \
$outer_seed \
$sample_out \
$initial_size \
$step_size \
$n_sim \
$stan_chains \
$stan_cores \
$stan_warmup \
$stan_iter \
$mapply_cores \
$pop_rho \
$model

end=$(date +%s)
total_time_secs=$(($end-$start))
c="${model}_${pop_rho}_${initial_size}_${sample_out}_${step_size}_${n_sim}_${stan_chains}_${stan_iter}"
echo $total_time_secs > $run_dir/output/elapsed_time_"$c".txt
#echo $total_time_secs > ~/output/elapsed_time_"$c".txt

# move output to project storage
#total_time_secs > $SNIC_TMP/elapsed_time_"$c".txt
# rsync -av ~/output /proj/snic2021-22-887/bayesian_correlations/results/
