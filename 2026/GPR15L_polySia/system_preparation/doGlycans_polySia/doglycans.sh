DGPATH=/home/courtade/doglycans
doglycans=/home/courtade/miniforge3/envs/doglycans/bin/prepreader

num=4

cat<< EOF> sequence_file.seq
polySia=ROH-(O1,C2)-{8SA-(O8,C2,<a>)}${num}-0SA
a=O8,C2,[C9 C8 O8 C2 90 C8 O8 C2 C1 90]
EOF


$doglycans \
	-f top/glycam06h.itp \
	-p top/GLYCAM_06h.prep \
	-s sequence_file.seq \
	-c POL_${num}8SA_0SA.pdb \
	-o POL_${num}8SA_0SA.itp \
	 --amber


