#!/bin/bash
sbatch <<EOT
#!/bin/sh

#SBATCH -A snic2020-15-75
#SBATCH -J $1 # name
#SBATCH -o $1_%j.out
#SBATCH -t $2 # time
#SBATCH -e $1\_%j.out
#SBATCH -p $3 # core or devel
#SBATCH -n $4 # cores
#SBATCH --mail-user=lecka48@liu.se
#SBATCH --mail-type=ALL

module purge
module load R/4.0.0
module --ignore-cache load gcc/10.1.0
module load openmpi/3.1.6
mpirun -n $4 R --no-save --args $6 $7 < $5 # 5 is scriptname
EOT