#!/bin/sh
# 2019.03.19
# Yong

source /lindorffgrp-isilon/wyong/software/GMX514/bin/GMXRC
gmx=/lindorffgrp-isilon/wyong/software/GMX514/bin/gmx_mpi_plumed241

theta=1.10

# minimization
$gmx grompp -f minimization.mdp -p all_LPMO_theta$theta.top -c LPMO_SOL_IONS.gro -o min.tpr -maxwarn 2
$gmx mdrun -deffnm min -v -ntomp 1

# relaxation
$gmx grompp -f relax.mdp -p all_LPMO_theta$theta.top -c min.gro -o md1.tpr -maxwarn 2
$gmx mdrun -deffnm md1 -v -ntomp 8

# productive
$gmx grompp -f md.mdp -p all_LPMO_theta$theta.top -c md1.gro -o md2.tpr -maxwarn 2
ln -s $gmx M3_LPMO_theta$theta
nohup ./M3_LPMO_theta$theta mdrun -deffnm md2 -v -ntomp 8 &

