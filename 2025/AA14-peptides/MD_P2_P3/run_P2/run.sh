#!/bin/bash
#
#SBATCH --job-name=pep-struc
#SBATCH --mem=2G
#SBATCH --cpus-per-gpu=4
#SBATCH --gpus=1

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/triumvirate/apps/plumed-2.9.0/lib
#export PLUMED_KERNEL=/triumvirate/apps/plumed-2.9.0/lib/libplumedKernel.so
#plumed=/triumvirate/apps/plumed-2.9.0/bin/plumed
#gmx=/triumvirate/apps/gromacs-2023_CUDA-12.0.0_plumed-2.9.0/bin/gmx
#gmx=gmx
gmx=/triumvirate/apps/gromacs2022.3_plumed2.8.3_cuda_mpi/bin/gmx_mpi

rm \#*
rm *tpr

STARTPDB=model1_NoCys.pdb
FFDIR=a99SBdisp.ff
FF=a99SBdisp
watergro=a99SBdisp.ff/a99SBdisp_water.gro


#generate topology
echo 1 | $gmx pdb2gmx -f $STARTPDB -o processed.gro -ignh -ff $FF -water select

# AMBER99SB-ILDN
#echo 7 | $gmx pdb2gmx -f $STARTPDB -o processed.gro -ignh -water tip3p


#solvation - define box + add water
$gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt dodecahedron
$gmx solvate -cp newbox.gro -cs $watergro -o solv.gro -p topol.top

# set up ions
$gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo 13 | $gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

#energy minimization
$gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
$gmx mdrun -v -deffnm em 

# generate energy file
echo 11 0 | $gmx energy -f em.edr -o potential.xvg  

# NVT Equilibration (stabilize temperature)
$gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
$gmx mdrun -v -deffnm nvt

#generate file for temperature progression
echo 17 0 |$gmx energy -f nvt.edr -o temperature.xvg

# Equilibration under NPT ensemble  (to stabilize pressure and density)
$gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
$gmx mdrun -v -deffnm npt

#generate file for pressure and density progression
echo 19 0 |$gmx energy -f npt.edr -o pressure.xvg
echo 25 0 |$gmx energy -f npt.edr -o density.xvg

# Run production MD for data collection
$gmx grompp -f md.mdp -c  npt.gro -t npt.cpt -p topol.top -o md.tpr
$gmx mdrun -v -deffnm md > gmx.log 2>&1
