
rm \#*
rm *tpr

STARTPDB=GPR15LG.pdb
FF=amber99sb-ildn
WATER=tip3p

gmx=gmx


#generate topology
$gmx pdb2gmx -f $STARTPDB -o processed.gro -ignh -ff $FF -water $WATER -ss
$gmx select -s processed.gro -select "resname CYS and name SG" -oi CYS.dat -on CYS.ndx
$gmx genrestr -f processed.gro -n CYS.ndx -disre -cutoff 0.4
# remember to add posres here

$gmx grompp -f mdp/minim.mdp -c pro_wSS.gro -p topol.top -o minim.tpr -maxwarn 2
$gmx mdrun -v -deffnm minim -nt 4


# remember to remove header and footer of generated topol.top and rename to prot.itp
