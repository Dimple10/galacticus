#!/bin/bash
#PBS -N MDPL_mcmc
#PBS -l nodes=4:ppn=16
#PBS -j oe
#PBS -o MDPL_mcmc2.log
#SBATCH --mem-per-cpu=8G
#PBS -m ea
#PBS -M dsarnaaik@carnegiescience.edu
#PBS -V
if [ ! -z ${PBS_O_WORKDIR+x} ]; then
 cd $PBS_O_WORKDIR
elif [ ! -z ${SLURM_SUBMIT_DIR+x} ]; then
 cd $SLURM_SUBMIT_DIR
fi

export LD_LIBRARY_PATH=/home/abenson/Galacticus/Tools/lib:/home/abenson/Galacticus/Tools/lib64:$LD_LIBRARY_PATH
export PATH=/home/abenson/Galacticus/Tools/bin:$PATH
export GALACTICUS_FCFLAGS="-fintrinsic-modules-path /home/abenson/Galacticus/Tools/finclude -fintrinsic-modules-path /home/abenson/Galacticus/Tools/include -fintrinsic-modules-path /home/abenson/Galacticus/Tools/include/gfortran -fintrinsic-modules-path /home/abenson/Galacticus/Tools/lib/gfortran/modules -L/home/abenson/Galacticus/Tools/lib -L/home/abenson/Galacticus/Tools/lib64 -L/home/abenson/Galacticus/Tools/lib -L/home/abenson/Galacticus/Tools/lib64"
export GALACTICUS_CFLAGS="-I/home/abenson/Galacticus/Tools/include"
export GALACTICUS_CPPFLAGS="-I/home/abenson/Galacticus/Tools/include"
ulimit -t unlimited
ulimit -c unlimited
export OMP_NUM_THREADS=1
/usr/bin/time -v mpirun --n 64 --map-by node --bind-to none -mca btl ^openib ./Galacticus.exe constraints/pipelines/darkMatter/progenitorMassFunctionConfig_MDPLTest2.xml
exit
