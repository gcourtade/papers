#!/bin/bash

py37="/groups/sbinlab/courtade/.conda/envs/py37/bin/python" 
export PATH=$PATH:/groups/sbinlab/software/openmpi-2.1.0/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/lustre/hpc/sbinlab/software/miniconda3/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/groups/sbinlab/software/openmpi-2.1.0/lib
gmx=/groups/sbinlab/wyong/usr/local/GMX2018.2/bin/gmx_mpi_plumed
export GMX_ALLOW_CPT_MISMATCH=1
export PLUMED_USE_LEPTON=yes
export PLUMED_KERNEL="/groups/sbinlab/wyong/usr/local/PLUMED250dev_PathCV/lib/libplumedKernel.so"
export LD_LIBRARY_PATH="/groups/sbinlab/wyong/usr/local/PLUMED250dev_PathCV/lib":$LD_LIBRARY_PATH



#echo 1 | $gmx gyrate -f md_prot.xtc -s md2.tpr -o rg.xvg
#$py37 ~/scripts/python3-block_analysis/average_block.py rg.xvg 1 | tee rg.stat

$gmx distance -f md_prot.xtc -s md2.tpr -select 10 -n index.ndx -oav dist.xvg
#$py37 ~/scripts/python3-block_analysis/average_block.py dist.xvg 1 | tee dist.stat

