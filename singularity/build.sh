#!/bin/bash

cd $PWD/build_dir
sudo singularity build bayesian_correlations.sif bayesian_correlations.def
cd ..
