#!/bin/bash
#
#SBATCH --job-name=AA14_dimer
#SBATCH --mem=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=10-00:00:00
#SBATCH --gpus=1

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/triumvirate/apps/plumed-2.9.0/lib
#export PLUMED_KERNEL=/triumvirate/apps/plumed-2.9.0/lib/libplumedKernel.so
#plumed=/triumvirate/apps/plumed-2.9.0/bin/plumed
#gmx=/triumvirate/apps/gromacs-2023_CUDA-12.0.0_plumed-2.9.0/bin/gmx
gmx=/triumvirate/apps/gromacs2022.3_plumed2.8.3_cuda_mpi/bin/gmx_mpi
rm \#*
rm *tpr

STARTPDB=af3model0.pdb
FFDIR=a99SBdisp-CU.ff
FF=a99SBdisp-CU
watergro=a99SBdisp-CU.ff/a99SBdisp_water.gro

cp start.top topol.top
#solvation - define box + add water
$gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt dodecahedron
$gmx solvate -cp newbox.gro -cs $watergro -o solv.gro -p topol.top

# set up ions
$gmx grompp -f mdp/ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo 13 | $gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1


#energy minimization
$gmx grompp -f mdp/emin.mdp -c solv_ions.gro -p topol.top -o em.tpr
$gmx mdrun -v -deffnm em


# generate energy file
#echo 11 0 | $gmx energy -f em.edr -o potential.xvg  

# NVT Equilibration (stabilize temperature)
$gmx grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
$gmx mdrun -v -deffnm nvt

#generate file for temperature progression
#echo 17 0 |$gmx energy -f nvt.edr -o temperature.xvg

# Equilibration under NPT ensemble  (to stabilize pressure and density)
$gmx grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
$gmx mdrun -v -deffnm npt

#generate file for pressure and density progression
#echo 19 0 |$gmx energy -f npt.edr -o pressure.xvg
#echo 25 0 |$gmx energy -f npt.edr -o density.xvg

# Run production MD for data collection
$gmx grompp -f mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
$gmx mdrun -v -deffnm md
