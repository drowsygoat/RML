#!/bin/bash

#SBATCH -A snic2020-15-75
#SBATCH -J 12_dec
#SBATCH -o 12_dec_%j.out
#SBATCH -t 01:00:00
#SBATCH -e 12_dec_%j.out
#SBATCH -p devel
#SBATCH -n 3
#SBATCH --mail-user=lecka48@liu.se
#SBATCH --mail-type=ALL

module purge
module load R/4.0.0
module --ignore-cache load gcc/10.1.0
module load openmpi/3.1.6
mpirun -n 3 R --no-save < 12_dec_diamond.R $1
