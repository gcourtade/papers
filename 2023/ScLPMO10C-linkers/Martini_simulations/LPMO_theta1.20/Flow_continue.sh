#!/bin/bash
#
#SBATCH --job-name=1.20-CelS2
#SBATCH --nodes=1                # number of nodes
#SBATCH --ntasks=8
#SBATCH --partition=sbinlab_ib
#SBATCH --time=24:00:00

export PATH=$PATH:/groups/sbinlab/software/openmpi-2.1.0/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/hpc/sbinlab/software/miniconda3/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/groups/sbinlab/software/openmpi-2.1.0/lib
gmx=/groups/sbinlab/wyong/usr/local/GMX2018.2/bin/gmx_mpi_plumed
export GMX_ALLOW_CPT_MISMATCH=1
export PLUMED_USE_LEPTON=yes
export PLUMED_KERNEL="/groups/sbinlab/wyong/usr/local/PLUMED250dev_PathCV/lib/libplumedKernel.so"
export LD_LIBRARY_PATH="/groups/sbinlab/wyong/usr/local/PLUMED250dev_PathCV/lib":$LD_LIBRARY_PATH

mpirun -n 8 $gmx mdrun -deffnm md2 -v -ntomp 1 -noappend -cpi



