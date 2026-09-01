
tleap=/home/courtade/miniforge3/envs/ambertools/bin/tleap
acpype=/home/courtade/miniforge3/envs/ambertools/bin/acpype


cat <<EOF > tleap.in
source leaprc.GLYCAM_06j-1
source leaprc.protein.ff14SB   
source leaprc.lipid21
 
glycoprot = loadpdb bilayer_amber_compatible_H.pdb

saveamberparm glycoprot glycoprot.prmtop glycoprot.inpcrd 
savepdb glycoprot glycoprot.pdb
quit
EOF

$tleap -f tleap.in
$acpype -p glycoprot.prmtop -x glycoprot.inpcrd 

