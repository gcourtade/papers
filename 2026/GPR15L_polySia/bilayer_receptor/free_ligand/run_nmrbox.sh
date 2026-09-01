#!/bin/bash
#SBATCH --job-name=GPR15LR
#SBATCH --mem=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus=1


# number of CPU cores used by PLUMED
APP_PATH=/home/nmrbox/gcourtade
LIBTORCH=$APP_PATH/libtorch
export CPATH=${LIBTORCH}/include/torch/csrc/api/include/:${LIBTORCH}/include/:${LIBTORCH}/include/torch:$CPATH
export INCLUDE=${LIBTORCH}/include/torch/csrc/api/include/:${LIBTORCH}/include/:${LIBTORCH}/include/torch:$INCLUDE
export LIBRARY_PATH=${LIBTORCH}/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=${LIBTORCH}/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$APP_PATH/plumed2-libtorch/lib
export PLUMED_KERNEL=$APP_PATH/plumed2-libtorch/lib/libplumedKernel.so
plumed=$APP_PATH/plumed2-libtorch/bin/plumed
gmx=$APP_PATH/gromacs2024.3_plumed2libtorch_cuda/bin/gmx

export CUDA_VISIBLE_DEVICES=0


#gmx grompp -f mdp/md.mdp -c npt.gro -p topol.top -pp all.top -o md.tpr
gmx mdrun -v -deffnm md   -ntomp 16 -cpi

