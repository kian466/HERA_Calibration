#!/lustre/aoc/projects/hera/kshahin/miniconda3/envs/hera3/bin/bash
#PBS -N Calibration
#PBS -q hera
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:10:00
#PBS -l vmem=16g
#PBS -j oe
#PBS -o Calibration.out
#PBS -m be
#PBS -M kian.shahin@berkeley.edu
#PBS -t 1-576%32
source deactivate
source activate hera3
python /home/kshahin/HERA_Calibration/calibration.py ${PBS_ARRAYID} /lustre/aoc/projects/hera/aewallwi/time_dep_flagging/2458098 /lustre/aoc/projects/hera/aewallwi/tim\
e_dep_flagging/2458098

bash copy_output_to_lustre.sh


