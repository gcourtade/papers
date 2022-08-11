#!/bin/bash
#
#SBATCH --job-name=extended_backmap
#SBATCH --nodes=1                # number of nodes
#SBATCH --ntasks=1
#SBATCH --partition=sbinlab
#SBATCH --time=24:00:00

#echo wget http://cgmartini.nl/images/tools/backward/backward-v5.zip

xtc=md_prot.xtc
tpr=md2.tpr
AApdb=extended.pdb


gmx=/groups/sbinlab/software/gromacs-5.1.4/bin/gmx_mpi_plumed2.3.0
backmappingdir=/groups/sbinlab/courtade/software/backward-v5
FF=charmm27

# Step1: extract pdb frames
step=1
if [ $step -eq 1 ]; then
echo 1 | $gmx trjconv -f $xtc -o frame.pdb -pbc whole -skip 1 -sep -s $tpr 
fi

# Step2: prepare all-atom topology

if [ $step -eq 1 ]; then
$gmx pdb2gmx -f $AApdb -o AA.gro -ignh -water none -ff $FF
fi

rm -f \#*

# Step3: backmapping

#i=0
# starting frame = number of current frames
i=`ls AA*pdb | wc -l`
nframe=`ls -d frame*pdb | wc -l`

while [ $i -lt $nframe ]
do

rm -f \#*
rm backmapped.gro

name=frame${i}
CG=${name}.pdb

#echo $i $name
#echo bash $backmappingdir/initram-v5.sh -f $CGgro -p topol.top -to charmm36
#bash $backmappingdir/initram-v5.sh -f $CGgro -p topol.top -to charmm36 > /dev/null

# the quick version
bash $backmappingdir/initram-v5_150min.sh -f $CG -p topol.top -to charmm36

$gmx editconf -f backmapped.gro -o AA_${name}.pdb 
echo $gmx editconf -f backmapped.gro -o AA_${name}.pdb 

i=$(( $i + 1 ))

done
