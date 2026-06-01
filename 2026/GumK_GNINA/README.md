# Conformation-Dependent Donor Selectivity in the Xanthan Gum Glycosyltransferase GumK Revealed by AI-Based Ensemble Docking 


Authors: Davide Luciano, Knall Anna Tova Alenfalk, and Gaston Courtade1

To run the docking method use the following script in bash:

```bash
for i in {0..30}; do

mkdir rep_$i
cd rep_$i

gnina -r ../receptor.pdb -l ../ligand.mol2 --autobox_ligand ../receptor_autobox.pdb --exhaustiveness 16 --autobox_extend 1 --min_rmsd_filter 1 --num_modes 20 -o WT1_glcA.mol2 --seed $RANDOM > gnina.log 2>&1

cd ..

done
```

receptor is named WTB.pdb or WT1.pdb in the directoryes for confomration 0 or 1. ligand.mol2 is named depending on the specifci ligand. 
