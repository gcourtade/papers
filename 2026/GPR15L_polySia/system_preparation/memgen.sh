#!/bin/bash

memgen=/home/courtade/miniforge3/envs/ambertools/bin/packmol-memgen
python=/home/courtade/miniforge3/envs/ambertools/bin/python


cp start_files/ProtH.pdb Prot.pdb

$memgen --pdb Prot.pdb --lipids POPC --distxy_fix 100 --parametrize
$python -c "import parmed as pmd; parm = pmd.load_file('bilayer_Prot_lipid.top', 'bilayer_Prot_lipid.crd'); parm.save('bilayer_Prot_lipid_GMX.top', format='gromacs'); parm.save('bilayer_Prot_lipid_GMX.gro')"

