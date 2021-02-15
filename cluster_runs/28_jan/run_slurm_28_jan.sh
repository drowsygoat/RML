#!/bin/bash
sbatch <<EOT
#!/bin/sh

#SBATCH -A snic2020-15-75
#SBATCH -J 12_dec_dist
#SBATCH -o 12_dec_dist_%j.out
#SBATCH -t 01:00:00
#SBATCH -e 12_dec_dist_%j.out
#SBATCH -p devel
#SBATCH -n 10
#SBATCH --mail-user=lecka48@liu.se
#SBATCH --mail-type=ALL

module purge
module load R/4.0.0
module --ignore-cache load gcc/10.1.0
module load openmpi/3.1.6
mpirun -n 10 R --no-save --args $1 $2 < 12_dec_distances.R
EOT