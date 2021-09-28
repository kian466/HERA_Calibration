#!/bin/bash -l
#PBS -N Calibration
#PBS -q hera
#PBS -l nodes=1:ppn=1
#PBS -l walltime=03:00:00
#PBS -l vmem=16g
#PBS -j oe
#PBS -o Calibration_nonredundant.out
#PBS -m be
#PBS -M kian.shahin@berkeley.edu
#PBS -t 0-17

source ~/.bashrc
conda deactivate
conda activate hera3

python /users/kshahin/kshahin/HERA_Calibration/calibration_nonredundant_imperfect.py ${PBS_ARRAYID} /lustre/aoc/projects/hera/schoudhu/vis_redun_60time_90mjy.uvh5 /lustre/aoc/projects/hera/schoudhu/vis_ellip_60time.uvh5
