#!/bin/sh
# 2019.03.19
# Yong

export PATH=$PATH:/groups/sbinlab/software/openmpi-2.1.0/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/hpc/sbinlab/software/miniconda3/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/groups/sbinlab/software/openmpi-2.1.0/lib
gmx=/groups/sbinlab/wyong/usr/local/GMX2018.2/bin/gmx_mpi_plumed
export GMX_ALLOW_CPT_MISMATCH=1
export PLUMED_USE_LEPTON=yes
export PLUMED_KERNEL="/groups/sbinlab/wyong/usr/local/PLUMED250dev_PathCV/lib/libplumedKernel.so"
export LD_LIBRARY_PATH="/groups/sbinlab/wyong/usr/local/PLUMED250dev_PathCV/lib":$LD_LIBRARY_PATH

# install martinize by following https://github.com/marrink-lab/vermouth-martinize/blob/master/README.md
# my version is here:
martinize=/groups/sbinlab/courtade/.conda/envs/py37/bin/martinize2
#martinize=/lindorffgrp-isilon/wyong/software/miniconda3/bin/martinize2
#wget http://cgmartini.nl/images/tools/insane/insane.py
insane=insane.py

FF=martini304

dssp=/groups/sbinlab/courtade/.conda/envs/py37/bin/mkdssp

pdb=polyS.pdb
#$gmx pdb2gmx -f $pdb -o clean-$pdb -water none -his
#$gmx editconf -f $pdb -o clean-$pdb
cp $pdb clean-$pdb
ffdir=/lustre/hpc/sbinlab/courtade/simulations/Martini30b417

#==============
# CG
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
$martinize -f clean-$pdb -o AA_topol.top -x AA_CG.pdb -dssp $dssp -ff $FF -ed -cys auto -elastic -scfix -ff-dir $ffdir/v.3.0.4.17/martini3-protein/force_fields/ -map-dir $ffdir/v.3.0.4.17/martini3-protein/mappings/
# for IDP, we should remove -dssp and -elastic options
#python $martinize -f $pdb -o AA_topol.top -x AA_CG.pdb -dssp $dssp -ff $FF -ed -cys auto -elastic -scfix -ff-dir $ffdir/v.3.0.4.17/martini3-protein/force_fields/ -map-dir $ffdir/v.3.0.4.17/martini3-protein/mappings/

#=============
# add solvent/ions
$gmx editconf -f AA_CG.pdb -o AA_CG.gro -bt cubic -d 1
python2.7 $insane -f AA_CG.gro -o AA_SOL_IONS.gro -pbc keep -salt 0.1 -sol W -center -p AA_topol_SOL_IONS.top 

# modify top
perl -pi -e's/#include "martini.itp"//g' AA_topol_SOL_IONS.top
perl -pi -e's/NA\+/NA/g' AA_topol_SOL_IONS.top
perl -pi -e's/CL-/CL/g' AA_topol_SOL_IONS.top
mv molecule_0.itp AA.itp
perl -pi -e's/molecule_0/Protein/g' AA.itp
cat <<EOF > others.top
#include "$ffdir/v.3.0.4.17/martini_v3.0.4.itp"
#include "AA.itp"
#include "$ffdir/v.3.0.4.17/martini_v3.0_phospholipids.itp"
#include "$ffdir/v.3.0.4.17/martini_v3.0_ions.itp"
#include "$ffdir/v.3.0.4.17/martini_v3.0_solvents.itp"
EOF
cat others.top AA_topol_SOL_IONS.top >a
mv a AA_topol_SOL_IONS.top

#!!! for the simulations of multi-domain proteins
#!!! remember to visulize the elastic network model and remove the linker-domain and domain-domain contacts

