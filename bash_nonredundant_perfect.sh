#!/bin/bash -l
#PBS -N Calibration_Perfect
#PBS -q hera
#PBS -l nodes=1:ppn=1
#PBS -l walltime=03:00:00
#PBS -l vmem=16g
#PBS -j oe
#PBS -o Calibration_nonredundant_perfect.out
#PBS -m be
#PBS -M kian.shahin@berkeley.edu
#PBS -t 0-17

source ~/.bashrc
conda deactivate
conda activate hera3

python /users/kshahin/kshahin/HERA_Calibration/calibration_nonredundant_perfect.py ${PBS_ARRAYID} /lustre/aoc/projects/hera/schoudhu/vis_ellip_60time.uvh5 /lustre/aoc/projects/hera/schoudhu/vis_ellip_60time.uvh5
