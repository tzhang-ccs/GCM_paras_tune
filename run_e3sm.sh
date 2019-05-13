#! /bin/sh
# Batch system directives
#SBATCH  --job-name=low
#SBATCH  --nodes=40
#SBATCH  --output=low_ne30_ne30_test2222.run.%j 
#SBATCH  --exclusive 
#SBATCH  --constraint=knl,quad,cache
#SBATCH  -p regular
#SBATCH  --time 10:15:00

export 'OMP_STACKSIZE'='128M' 
export 'OMP_PROC_BIND'='spread'
export 'OMP_NUM_THREADS'='4'
export 'OMP_PLACES'='threads'
srun  --label  -n 2560 -c 4  --cpu_bind=cores  /global/cscratch1/sd/zhangtao/acme_scratch/cori-knl/low_ne30/bld/e3sm.exe 
