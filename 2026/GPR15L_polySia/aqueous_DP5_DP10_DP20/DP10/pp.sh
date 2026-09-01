echo 27 | gmx trjconv -f npt.gro -s md.tpr -pbc whole -o md_complex.gro -n pp_index.ndx
echo 27 | gmx trjconv -f md.xtc -s md.tpr -pbc whole -o md_complex.xtc -n pp_index.ndx
