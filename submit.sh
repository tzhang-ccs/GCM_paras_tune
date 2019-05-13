#! /bin/sh
# Batch system directives
#SBATCH  --job-name=e3sm
#SBATCH  --nodes=40
#SBATCH  --output=output.%j 
#SBATCH  --error=error.%j
##SBATCH  --exclusive 
#SBATCH  --constraint=knl,quad,cache
##SBATCH  -p scavenger
#SBATCH  -p regular
#SBATCH  --time 14:00:00
##SBATCH --time-min 1:30:00

./downhill_simplex
