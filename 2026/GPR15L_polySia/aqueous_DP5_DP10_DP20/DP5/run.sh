#!/bin/bash
#
#SBATCH --job-name=complex_DP5_work
#SBATCH --mem=2G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=07-00:00:00
#SBATCH --gpus=1
#SBATCH --nodelist=herod

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/triumvirate/apps/plumed2.8.3/lib
export PLUMED_KERNEL=/triumvirate/apps/plumed2.8.3/lib/libplumedKernel.so
plumed=/triumvirate/apps/plumed2.8.3/bin/plumed
gmx=/triumvirate/apps/gromacs2022.3_plumed2.8.3_cuda_mpi/bin/gmx_mpi

rm rep*

cat <<EOF >plumed.dat
gA: GROUP ATOMS=1-256
gB: GROUP ATOMS=938-1120

comA: COM ATOMS=gA
comB: COM ATOMS=gB

d: DISTANCE ATOMS=comA,comB

uw: UPPER_WALLS ARG=d AT=7.3 KAPPA=10000

metad: METAD ...
  ARG=d
  SIGMA=0.005
  #HEIGHT=0.04184 
  HEIGHT=0.4184 
  PACE=45
  FILE=HILLS
...


PRINT STRIDE=500 ARG=* FILE=COLVAR


EOF
NREPS=13
BASE_SEED=2007
START=3
END=$((START + NREPS - 1))

# starting conformation
echo 0 | $gmx trjconv -f md.cpt -s md.tpr -o bound.gro

for i in $(seq "$START" "$END"); do
  repdir="rep_${i}"
  seed=$((BASE_SEED + i))
  echo "===> Start ${repdir}"
  mkdir -p "$repdir"
  cp mdp/nvt_seeded.mdp rep_${i}/nvt_seeded.mdp
  cd "$repdir"
  sed -i 's/GEN_SEED/'$seed'/' nvt_seeded.mdp
  # NVT Equilibration (stabilize temperature)
  $gmx grompp -f nvt_seeded.mdp -c ../bound.gro -r ../bound.gro -p ../topol.top -o nvt.tpr -n ../index.ndx
  $gmx mdrun -v -deffnm nvt

  # Equilibration under NPT ensemble  (to stabilize pressure and density)
  $gmx grompp -f ../mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p ../topol.top -o npt.tpr -n ../index.ndx
  $gmx mdrun -v -deffnm npt

  # Run production MD for data collection
  $gmx grompp -f ../mdp/md_short.mdp -c npt.gro -t npt.cpt -p ../topol.top -o md.tpr -n ../index.ndx
  $gmx mdrun -v -deffnm md -plumed ../plumed.dat
  cd ..
done


