#!/bin/bash
#
#SBATCH --job-name=complex_polySia_DP20
#SBATCH --mem=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=14-00:00:00
#SBATCH --gpus=1

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/triumvirate/apps/plumed-2.9.0/lib
#export PLUMED_KERNEL=/triumvirate/apps/plumed-2.9.0/lib/libplumedKernel.so
#plumed=/triumvirate/apps/plumed-2.9.0/bin/plumed
#gmx=/triumvirate/apps/gromacs-2023_CUDA-12.0.0_plumed-2.9.0/bin/gmx
gmx=/triumvirate/apps/gromacs2022.3_plumed2.8.3_cuda_mpi/bin/gmx_mpi

$gmx mdrun -v -deffnm md -cpi
