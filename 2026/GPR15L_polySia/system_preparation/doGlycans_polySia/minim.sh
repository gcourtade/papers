gmx=gmx

cat<< EOF> topol.top
; Include forcefield parameters
#include "amber99sb-ildn.ff/forcefield.itp"
#include "top/glycam06h.itp"
#include "POL.itp"

; Include water topology
#include "amber99sb-ildn.ff/tip3p.itp"

; Include ions topology
#include "amber99sb-ildn.ff/ions.itp"

[ system ]
; Name
polySia in water

[ molecules ]
; Compound        #mols
polySia     5
EOF

#solvation - define box + add water
$gmx editconf -f POL_5.pdb -o newbox.gro -c -d 1.0 -bt dodecahedron
$gmx solvate -cp newbox.gro  -o solv.gro -p topol.top


# set up ions
$gmx grompp -f mdp/ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo 6 | $gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral


#energy minimization
$gmx grompp -f mdp/minim.mdp -c solv_ions.gro -p topol.top -o em.tpr
$gmx mdrun -v -deffnm em -nt 4

