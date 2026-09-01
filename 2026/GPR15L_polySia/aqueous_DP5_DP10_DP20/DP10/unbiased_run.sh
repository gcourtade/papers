#!/bin/bash
#
#SBATCH --job-name=complex_polySia_DP10
#SBATCH --mem=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=07-00:00:00
#SBATCH --gpus=1

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/triumvirate/apps/plumed-2.9.0/lib
#export PLUMED_KERNEL=/triumvirate/apps/plumed-2.9.0/lib/libplumedKernel.so
#plumed=/triumvirate/apps/plumed-2.9.0/bin/plumed
#gmx=/triumvirate/apps/gromacs-2023_CUDA-12.0.0_plumed-2.9.0/bin/gmx
gmx=/triumvirate/apps/gromacs2022.3_plumed2.8.3_cuda_mpi/bin/gmx_mpi

#rm \#*
#rm *tpr
#
#STARTPDB=complex.pdb
#
#cp start.top topol.top
#
##solvation - define box + add water
#$gmx editconf -f $STARTPDB -o newbox.gro -c -d 1.0 -bt dodecahedron
#$gmx solvate -cp newbox.gro -cs $watergro -o solv.gro -p topol.top
#
## set up ions
#$gmx grompp -f mdp/ions.mdp -c solv.gro -p topol.top -o ions.tpr
#echo 17 | $gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
#
##energy minimization
#$gmx grompp -f mdp/emin.mdp -c solv_ions.gro -p topol.top -o em.tpr
#$gmx mdrun -v -deffnm em
#
## on the first run exit and make_ndx for the three sugar residues, then comment exit and rerun
##exit
#
## NVT Equilibration (stabilize temperature)
#$gmx grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx
#$gmx mdrun -v -deffnm nvt
#
## Equilibration under NPT ensemble  (to stabilize pressure and density)
#$gmx grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -n index.ndx
#$gmx mdrun -v -deffnm npt
#
## Run production MD for data collection
#$gmx grompp -f mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -n index.ndx
$gmx mdrun -v -deffnm md -cpi
