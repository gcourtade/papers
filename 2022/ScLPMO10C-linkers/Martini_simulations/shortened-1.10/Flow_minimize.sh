#!/bin/bash
#
#SBATCH --job-name=1.10-shortened
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

rm *tpr

theta=1.10
$gmx grompp -f minimization.mdp -p AA_topol_SOL_IONS.top -c AA_SOL_IONS.gro -o min.tpr -pp all_AA.top -maxwarn 2

perl /lustre/hpc/sbinlab/courtade/scripts/PWrescaling_Martini3b0417.pl -I all_AA.top -T $theta -O all_LPMO_theta$theta.top


# minimization
$gmx grompp -f minimization.mdp -p all_LPMO_theta$theta.top -c AA_SOL_IONS.gro -o min.tpr -maxwarn 2
mpirun -n 1 $gmx mdrun -deffnm min -v

# relaxation
$gmx grompp -f relax.mdp -p all_LPMO_theta$theta.top -c min.gro -o md1.tpr -maxwarn 2
mpirun -n 8 $gmx mdrun -deffnm md1 -v

# productive
$gmx grompp -f md.mdp -p all_LPMO_theta$theta.top -c md1.gro -o md2.tpr -maxwarn 2
#ln -s $gmx ${name}_$theta
mpirun -n 8 $gmx  mdrun -deffnm md2 -v

