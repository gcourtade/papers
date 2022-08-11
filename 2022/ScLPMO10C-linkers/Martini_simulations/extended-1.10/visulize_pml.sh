#!/bin/sh
#

#itpfile=ENM_current.itp
itpfile=ENM_corrected.itp
nline=`cat $itpfile | wc -l`
echo $nline
echo "load extended.pdb, AA" >vis.pml
echo "load AA_CG.pdb, CG" >>vis.pml
for i in `seq 1 1 $nline`
do
bi=`awk "NR==$i" $itpfile | awk '{print $1}'`
bj=`awk "NR==$i" $itpfile | awk '{print $2}'`
#echo "$i $bi $bj"
echo "$i $bi $bj" | awk '{printf("dist dist%d, CG and index %4d, CG and index %4d \n",$1,$2,$3)}' >>vis.pml
done
