#!/bin/bash

source /lindorffgrp-isilon/wyong/software/GMX514/bin/GMXRC
gmx=/lindorffgrp-isilon/wyong/software/GMX514/bin/gmx_mpi_plumed241

#ln -s $gmx theta1.10 
nohup ./theta1.10 mdrun -deffnm md2 -v -ntomp 8 -cpi -noappend&
